ROOT_DIR= $(shell pwd)

TARGETS = toolkits/update/ppsp_directed toolkits/update/viterbi_directed toolkits/update/ppwp_directed toolkits/update/ppnp_directed toolkits/update/bfs_directed toolkits/update/wcc_directed \
			toolkits/update/ppsp_undirected toolkits/update/viterbi_undirected toolkits/update/ppwp_undirected toolkits/update/ppnp_undirected toolkits/update/bfs_undirected toolkits/update/wcc_undirected \
			toolkits/no/ppsp_directed \
			toolkits/no/ppsp_undirected \
			toolkits/ub/ppsp_directed \
			toolkits/ub/ppsp_undirected \
			toolkits/lb/ppsp_directed \
			toolkits/lb/ppsp_undirected \
			toolkits/bs/ppsp_directed \
			toolkits/bs/ppsp_undirected \
			toolkits/all/ppsp_directed toolkits/all/viterbi_directed toolkits/all/ppwp_directed toolkits/all/ppnp_directed toolkits/all/bfs_directed toolkits/all/wcc_directed \
			toolkits/all/ppsp_undirected toolkits/all/viterbi_undirected toolkits/all/ppwp_undirected toolkits/all/ppnp_undirected toolkits/all/bfs_undirected toolkits/all/wcc_undirected \
			toolkits/pnp/ppsp_directed toolkits/pnp/viterbi_directed toolkits/pnp/ppwp_directed toolkits/pnp/ppnp_directed toolkits/pnp/bfs_directed toolkits/pnp/wcc_directed \
			toolkits/pnp/ppsp_undirected toolkits/pnp/viterbi_undirected toolkits/pnp/ppwp_undirected toolkits/pnp/ppnp_undirected toolkits/pnp/bfs_undirected toolkits/pnp/wcc_undirected \
			toolkits/tripoline/ppsp_directed toolkits/tripoline/viterbi_directed toolkits/tripoline/ppwp_directed toolkits/tripoline/ppnp_directed toolkits/tripoline/bfs_directed toolkits/tripoline/wcc_directed \
			toolkits/tripoline/ppsp_undirected toolkits/tripoline/viterbi_undirected toolkits/tripoline/ppwp_undirected toolkits/tripoline/ppnp_undirected toolkits/tripoline/bfs_undirected toolkits/tripoline/wcc_undirected \

#MACROS= -D INDEX_UPPER -D INDEX_LOWER -D BIDIRECTIONAL -D STATIC -D UPDATE -D BFS_OPT

MACROS_PLAIN = -D SGRAPH_CORE
MACROS_UPDATE = -D UPDATE -D STATIC -D SGRAPH_CORE
MACROS_NO = -D STATIC -D SGRAPH_CORE
MACROS_UB = -D INDEX_UPPER -D STATIC -D SGRAPH_CORE
MACROS_LB = -D INDEX_UPPER -D INDEX_LOWER -D STATIC -D SGRAPH_CORE
MACROS_BS = -D INDEX_UPPER -D BIDIRECTIONAL -D STATIC -D SGRAPH_CORE
MACROS_ALL = -D INDEX_UPPER -D INDEX_LOWER -D BIDIRECTIONAL -D STATIC -D SGRAPH_CORE
MACROS_ALL_BFS = -D INDEX_UPPER -D INDEX_LOWER -D BIDIRECTIONAL -D STATIC  -D SGRAPH_CORE -D BFS_OPT
MACROS_TRIPOLINE = -D STATIC -D TRIPOLINE_CORE
MACROS_PNP = -D STATIC -D PNP_CORE
MACROS_PNP_BFS = -D STATIC -D PNP_CORE -D BFS_OPT

MPICXX= mpicxx
CXXFLAGS= -Wall -O3 -std=c++14 -g -fopenmp -march=native -I $(ROOT_DIR) -I deps
HEADERS= $(shell find . -name '*.hpp')

all: $(TARGETS)

toolkits/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_PLAIN) -o $@ $<

toolkits/update/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_UPDATE) -o $@ $<

toolkits/no/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_NO) -o $@ $<

toolkits/ub/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_UB) -o $@ $<

toolkits/lb/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_LB) -o $@ $<

toolkits/bs/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_BS) -o $@ $<

toolkits/all/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_ALL) -o $@ $<

toolkits/all/bfs_directed: toolkits/bfs_directed.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_ALL_BFS) -o $@ $<

toolkits/all/bfs_undirected: toolkits/bfs_undirected.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_ALL_BFS) -o $@ $<
					
toolkits/tripoline/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_TRIPOLINE) -o $@ $<

toolkits/pnp/%: toolkits/%.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_PNP) -o $@ $<

toolkits/pnp/bfs_directed: toolkits/bfs_directed.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_PNP_BFS) -o $@ $<

toolkits/pnp/bfs_undirected: toolkits/bfs_undirected.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_PNP_BFS) -o $@ $<

toolkits/pnp/wcc_directed: toolkits/wcc_directed.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_PNP_BFS) -o $@ $<

toolkits/pnp/wcc_undirected: toolkits/wcc_undirected.cpp $(HEADERS)
					$(MPICXX) $(CXXFLAGS) $(MACROS_PNP_BFS) -o $@ $<

clean:
					rm -f $(TARGETS)
