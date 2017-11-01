# -*- coding: utf-8 -*-

from networkx.readwrite import json_graph
import json

def save_json_file(filename,graph):
    g = graph
    g_json = json_graph.node_link_data(g)
    json.dump(g_json,open(filename,'w'),indent=2)

def read_json_file(filename):
    with open(filename) as f:
        js_graph = json.load(f)
    return json_graph.node_link_graph(js_graph)