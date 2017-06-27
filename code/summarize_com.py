
# summarize_com.py
# Hojin Jung
# Last modified: 2017-06-27

import networkx as nx
from operator import itemgetter
from collections import defaultdict

def swap(a,b):
    if int(a) < int(b):
        return b,a
    return a,b

# summarize each community 
class SUM_COM:
    
    def __init__(self, e2c):
        self.e2c = e2c
        self.cid2edges, self.cid2nodes = defaultdict(set), defaultdict(set)
        self.edge2cid = {}
        self.initialize_c2e_c2n()

    # initialize 'community to edges (cid2edge)' and 'community to nodes (cid2nodes)'
    def initialize_c2e_c2n(self):
        for edge, cid in self.e2c.iteritems():
            self.cid2edges[cid].add(edge)
            self.cid2nodes[cid] |= set(edge)
        self.cid2edges, self.cid2nodes = dict(self.cid2edges), dict(self.cid2nodes)
    
    # summarization
    def summarize(self):
        for cid, edges in self.cid2edges.iteritems():
            
            numofnodes = len(self.cid2nodes[cid])
            numofedges = len(self.cid2edges[cid])

            # community: tree or star
            if numofedges == (numofnodes - 1):
                pass
            
            # community: clique
            elif numofedges == (((numofnodes - 1) * numofnodes) / 2):
                pass

            # community: other cases
            else:
                self.cluster_coef(cid, edges)
            
            # edge to community_id: (edge2cid)
            for e in self.cid2edges[cid]:
                # when edge exists
                if len(self.cid2nodes[cid]) == 1:
                    pass
                else:
                    e = swap(e[0], e[1]) # if neccessary
                    self.edge2cid[e] = cid
            
        return self.edge2cid

    # summarize other communities not tree, star, and clique
    def cluster_coef(self, cid, edges):
        
        G = nx.parse_edgelist(["%s %s" %(j[0], j[1]) for i, j in enumerate(edges)])

        # RuntimeWarning: use 'python -W ignore'
        # calculate average local clustering coefficient of a community
        avg_cluster = nx.average_clustering(G)
        
        # in case that you want to check a local clustering value of each node
        # cluster_list = nx.clustering(G)

        round_cluster = round(avg_cluster)
        
        t = 5 # Top_T
        
        if round_cluster == 0:
            Revised_G = self.remove_node_zero(G, top_t = t)

        # when avg_cluster is close to 1: density(high), high possibility to be assortative
        elif round_cluster == 1:
            Revised_G = self.remove_node_one(G, top_t = t)

        del self.cid2edges[cid]
        self.cid2edges[cid] = Revised_G.edges()
        del self.cid2nodes[cid]
        self.cid2nodes[cid] = Revised_G.nodes()
        

    def remove_node_zero(self, G, top_t = 0):
        # sort node in reverse degree order
        # ** Finding) Hubs are connected through nodes having as many degrees as the number of hubs(top_t)
        temp_t = top_t
        for node, deg in sorted(G.degree().items(), key = itemgetter(1), reverse = True):
            temp_t -= 1
            if temp_t <= 0 and deg != top_t:
                G.remove_node(node)
        return G
    
    def remove_node_one(self, G, top_t = 0):
        # nodes with high local clustering coefficients are alive (top_t)
        for node, cluster in sorted(nx.clustering(G).items(), key = itemgetter(1), reverse = True):
            top_t -= 1
            if top_t <= 0:
                G.remove_node(node)
        return G
