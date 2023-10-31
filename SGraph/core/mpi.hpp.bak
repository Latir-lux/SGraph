#ifndef MPI_HPP
#define MPI_HPP

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

template <typename T>
MPI_Datatype get_mpi_data_type()
{
	if (std::is_same<T, char>::value)
	{
		return MPI_CHAR;
	}
	else if (std::is_same<T, bool>::value)
	{
		return MPI_CHAR;
	}
	else if (std::is_same<T, unsigned char>::value)
	{
		return MPI_UNSIGNED_CHAR;
	}
	else if (std::is_same<T, int>::value)
	{
		return MPI_INT;
	}
	else if (std::is_same<T, unsigned>::value)
	{
		return MPI_UNSIGNED;
	}
	else if (std::is_same<T, long>::value)
	{
		return MPI_LONG;
	}
	else if (std::is_same<T, unsigned long>::value)
	{
		return MPI_UNSIGNED_LONG;
	}
	else if (std::is_same<T, float>::value)
	{
		return MPI_FLOAT;
	}
	else if (std::is_same<T, double>::value)
	{
		return MPI_DOUBLE;
	}
	else
	{
		printf("type not supported\n");
		exit(-1);
	}
}

class MPI_Instance
{
public:
	MPI_Comm index_comm, query_comm;
	MPI_Instance(int *argc, char ***argv)
	{
		int provided;
		MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
#ifdef PRINT_DEBUG_MESSAGES
		if (partition_id == 0)
		{
			printf("thread support level provided by MPI: ");
			switch (provided)
			{
			case MPI_THREAD_MULTIPLE:
				printf("MPI_THREAD_MULTIPLE\n");
				break;
			case MPI_THREAD_SERIALIZED:
				printf("MPI_THREAD_SERIALIZED\n");
				break;
			case MPI_THREAD_FUNNELED:
				printf("MPI_THREAD_FUNNELED\n");
				break;
			case MPI_THREAD_SINGLE:
				printf("MPI_THREAD_SINGLE\n");
				break;
			default:
				assert(false);
			}
		}
#endif
		MPI_Group world_group;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group);
		MPI_Comm_create(MPI_COMM_WORLD, world_group, &index_comm);
		MPI_Comm_create(MPI_COMM_WORLD, world_group, &query_comm);
	}
	~MPI_Instance()
	{
		MPI_Finalize();
	}
};

#endif
