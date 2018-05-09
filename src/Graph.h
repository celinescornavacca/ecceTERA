#ifndef GRAPH_H_
#define GRAPH_H_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };

namespace boost {
        BOOST_INSTALL_PROPERTY(vertex, properties);
        BOOST_INSTALL_PROPERTY(edge, properties);
}


template < typename VERTEXPROPERTIES, typename EDGEPROPERTIES >
class Graph
{

	protected:

        /* an adjacency_list like we need it */
        typedef boost::adjacency_list<
                boost::vecS, 
                boost::vecS, // vertex container
                boost::bidirectionalS, // directed graph
                boost::property<vertex_properties_t, VERTEXPROPERTIES>,
                boost::property<edge_properties_t, EDGEPROPERTIES>
        > GraphContainer;

	
		GraphContainer graph;
		
	public:

        /* a bunch of graph-specific typedefs */
        typedef typename boost::graph_traits<GraphContainer>::vertex_descriptor
            Vertex;
        typedef typename boost::graph_traits<GraphContainer>::edge_descriptor 
            Edge;
        typedef typename boost::graph_traits<GraphContainer>::vertex_iterator 
            vertex_iter;
        typedef typename boost::graph_traits<GraphContainer>::edge_iterator 
            edge_iter;
        typedef typename boost::graph_traits<GraphContainer>::adjacency_iterator
            adjacency_iter; 
        typedef typename boost::graph_traits<GraphContainer>::in_edge_iterator
            in_edge_iter;


        typedef std::pair<adjacency_iter, adjacency_iter> 
            adjacency_vertex_range_t;
    	typedef std::pair<in_edge_iter, in_edge_iter> in_edge_range_t;

        typedef std::pair<vertex_iter, vertex_iter> vertex_range_t;
        typedef std::pair<edge_iter, edge_iter> edge_range_t;


        /* constructors etc. */
        Graph() {}
        Graph(const Graph& g) : graph(g.graph) {}
        virtual ~Graph() {}


        /* structure modification methods */
        Vertex addVertex(const VERTEXPROPERTIES& prop) {
           	Vertex v = add_vertex(graph);
        	properties(v) = prop;
            return v;
        }
        void removeVertex( const Vertex& v ) {
            clear_vertex(v, graph);
            remove_vertex(v, graph);
        }

		bool isEdgeInGraph(const Vertex& v1, const Vertex& v2) {
			return edge(v1, v2, graph).second;

		}
        
         void addEdge(const Vertex& v1, const Vertex& v2, 
                 const EDGEPROPERTIES& prop_12) {
			if(!isEdgeInGraph(v1,v2)){
				Edge addedEdge1 = add_edge(v1, v2, graph).first;
              	properties(addedEdge1) = prop_12;
            }            
        }
	

        /* property access */
        VERTEXPROPERTIES& properties(const Vertex& v) {
           return (get(vertex_properties, graph))[v];
        }
        const VERTEXPROPERTIES& properties(const Vertex& v) const {
            return (get(vertex_properties, graph))[v];
        }
        EDGEPROPERTIES& properties(const Edge& e) {
            return (get(edge_properties, graph))[e];
        }
        const EDGEPROPERTIES& properties(const Edge& e) const {
            return (get(edge_properties, graph))[e];
        }


        /* selectors and properties */
        const GraphContainer& getGraph() const { return graph; }
		edge_range_t getEdges() const { return edges(graph); }
        
        VERTEXPROPERTIES getEdgeSource(const Edge& e) {
 			return properties(source(e, graph));
        }
        VERTEXPROPERTIES getEdgeTarget(const Edge& e) {
 			return properties(target(e, graph));
        }	  		
  		
		Vertex returnSource(const Edge& v) { return source(v, graph); }
        vertex_range_t getVertices() const { return vertices(graph); }

        adjacency_vertex_range_t getAdjacentVertices(const Vertex& v) const
        {
        	return adjacent_vertices(v, graph);
        }
        
        in_edge_range_t getInEdges(const Vertex& v) const {
        	return in_edges(v, graph);
        }
        int getVertexCount() const { return num_vertices(graph); }
        int getVertexOutDegree(const Vertex& v) const {
            return out_degree(v, graph);
        }
        int getVertexInDegree(const Vertex& v) const {
            return in_degree(v, graph);
        }



        /* operators */
        Graph& operator=(const Graph &rhs)
        {
            graph = rhs.graph;
            return *this;
        }

};

#endif
