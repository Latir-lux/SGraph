#ifndef SCHEDULER_HPP
#define SCHEDULER_HPP

#include "core/system.hpp"
#include <typeinfo>

template <typename Weight, typename EdgeData, typename Init, typename Forward, typename Backward, typename Merge>
class Scheduler
{
public:
    int argc;
    char **argv;
    int query_threads = 48;
    int index_threads = 48;

    bool symmetric;
    VertexId vertices;
    std::vector<EdgeUnit<EdgeData>> init_edges;

    uint64_t active_num = 0;
    uint64_t iteration_num = 0;
    uint64_t active_on_iteration[10000] = {0};
    int compute_num = 0;
    double compute_time = 0;

    Scheduler(int _argc, char **_argv, bool _symmetric)
    {
        argc = _argc;
        argv = _argv;
        symmetric = _symmetric;
    }
    VertexId get_random_vertex(VertexId vertices)
    {
        return 1.0 * rand() / RAND_MAX * vertices;//这句话会返回一个0~vertices之间的随机数
    }
    void compute(System<Weight, EdgeData, Init, Forward, Backward, Merge> &system)
    {
        VertexId source = get_random_vertex(vertices);
        VertexId sink = get_random_vertex(vertices);//终点
        uint64_t cur_active_num = 0;
        uint64_t cur_iteration_num = 0;
#ifdef PNP_CORE
        Weight res = system.pnp(source, sink, cur_active_num, cur_iteration_num, active_on_iteration);
#endif
#ifdef TRIPOLINE_CORE
        Weight res = system.tripoline(source, sink, cur_active_num, cur_iteration_num, active_on_iteration);
#endif
#ifdef SGRAPH_CORE
        Weight res = system.compute(source, sink, cur_active_num, cur_iteration_num, active_on_iteration);
#endif
        if (system.graph.partition_id == 0)
        {
            if (typeid(Weight) == typeid(float))
                printf("source = %d sink = %d res = %f\n", source, sink, res);
            else if (typeid(Weight) == typeid(int))
                printf("source = %d sink = %d res = %d\n", source, sink, res);
            else if (typeid(Weight) == typeid(bool))
                printf("source = %d sink = %d res = %d\n", source, sink, res);
            printf("active num = %lu iteration num = %lu\n", cur_active_num, cur_iteration_num);
        }
        active_num += cur_active_num;
        iteration_num += cur_iteration_num;
    }

    void work_static()
    {
        srand(0);
        MPI_Instance mpi(&argc, &argv);
        if (argc < 3)
        {
            printf("ppsp [file] [vertices]\n");
            exit(-1);
        }
        std::string file = argv[1];
        vertices = std::atoi(argv[2]);

        System<Weight, EdgeData, Init, Forward, Backward, Merge> system(mpi, vertices, symmetric, 48, 48);//后面两个参数是query_threads和index_threads
        omp_set_num_threads(96);

        system.load_file(file); // In this version of code, we use some tricks to reduce memory cost so that we can
                                // evaluate larger graphs with limited memory. But this is only for evaluation. The
                                // decoupled system architecture doesn't work. To let it work, disable compile macro
                                // `STATIC` and mutate graph like:
        // {
        //     std::vector<EdgeUnit<EdgeData>> edges;
        //     get_edge_vector(file, edges);
        //     system.add_edges(edges);
        //     system.step();
        // }

        omp_set_num_threads(48);

#ifndef PNP_CORE
        MPI_Barrier(MPI_COMM_WORLD);
        double t = -get_time();//一个巧妙的方法，用一个变量完成时间统计
        system.build_index();
        MPI_Barrier(MPI_COMM_WORLD);
        t += get_time();
        if (system.graph.partition_id == 0)
            printf("build time = %f\n", t);
        system.step();//step函数会将更新写入Adjlist中
#endif

        for (int i = 0; i < 500; i++)
        {
            MPI_Barrier(mpi.query_comm);
            double t = -get_time();
            compute(system);
            MPI_Barrier(mpi.query_comm);
            t += get_time();
            if (system.graph.partition_id == 0)
                printf("compute time = %f\n", t);
            compute_time += t;
            compute_num++;
        }

        if (system.graph.partition_id == 0)
        {
            printf("\n");
            printf("compute_num = %d\n", compute_num);
            printf("compute_time = %f\n", compute_time);
            printf("avg time = %f\n", compute_time / compute_num);
            printf("\n");
        }
        if (system.graph.partition_id == 0)
        {
            printf("compute_num = %d\n", compute_num);
            printf("active_num = %lu\n", active_num);
            printf("avg active = %f\n", (double)active_num / compute_num);
            printf("\n");
        }
        if (system.graph.partition_id == 0)
        {
            printf("compute_num = %d\n", compute_num);
            printf("iteration_num = %lu\n", iteration_num);
            printf("avg iteration = %f\n", (double)iteration_num / compute_num);
            printf("\n");
        }
        if (system.graph.partition_id == 0)
        {
            for (int i = 0; i < 10000; i++)
            {
                if (active_on_iteration[i] == 0)
                    break;
                printf("iteration %d active %f\n", i, (double)active_on_iteration[i] / compute_num);//active_on_iteration统计每一轮活跃顶点数，为了减小统计误差，执行了compute_num次，所以最终结果要除以compute_num。
            }
        }
    }

    void work_update()
    {//updata只做图更新工作,并不进行计算
        srand(0);
        MPI_Instance mpi(&argc, &argv);
        if (argc < 3)
        {
            printf("ppsp [file] [vertices]\n");
            exit(-1);
        }
        std::string file = argv[1];
        vertices = std::atoi(argv[2]);

        System<Weight, EdgeData, Init, Forward, Backward, Merge> system(mpi, vertices, symmetric, 48, 48);//后面两个参数是query_threads和index_threads
        omp_set_num_threads(48);

        EdgeUnit<EdgeData> **dynamic_edges_add;
        EdgeUnit<EdgeData> **dynamic_edges_del;
        system.load_file_update(file, dynamic_edges_add, dynamic_edges_del);

        MPI_Barrier(MPI_COMM_WORLD);
        double t = -get_time();
        system.build_index();
        MPI_Barrier(MPI_COMM_WORLD);
        t += get_time();
        if (system.graph.partition_id == 0)
            printf("build time = %f\n", t);
        system.step();

        double tot_t[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++)
        {
            int array_size;
            if (i == 0)
                array_size = 5000;
            if (i == 1)
                array_size = 50000;
            if (i == 2)
                array_size = 500000;
            for (int j = 0; j < 100; j++)
            {
                system.sim_update(dynamic_edges_add[i * 100 + j], dynamic_edges_del[i * 100 + j], array_size);
                MPI_Barrier(mpi.index_comm);
                double t = -get_time();
                system.modify_index();
                MPI_Barrier(mpi.index_comm);
                t += get_time();
                tot_t[i] += t;
                if (system.graph.partition_id == 0)
                    printf("array_size = %d modify_index time = %f\n", array_size, t);
            }
        }

        if (system.graph.partition_id == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                int array_size;
                if (i == 0)
                    array_size = 5000;
                if (i == 1)
                    array_size = 50000;
                if (i == 2)
                    array_size = 500000;
                printf("update %d edges per second\n", array_size * 2);
                printf("modify_index_num = %d\n", 100);
                printf("modify_index_time = %f\n", tot_t[i]);
                printf("avg time = %f\n", tot_t[i] / 100);
            }
        }
    }

};

#endif
