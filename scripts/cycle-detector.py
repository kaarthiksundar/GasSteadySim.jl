import networkx as nx
import json


def read_data(folder: str):
    file = folder + 'network.json'
    f = open(file)
    data = json.load(f)
    f.close()
    return data

def create_graph(data: dict):
    G = nx.Graph()
    G.add_nodes_from([int(i) for i in data["nodes"]])
    G.add_edges_from([(int(val["from_node"]), int(val["to_node"])) for key, val in data["compressors"].items()])
    # G.add_edges_from([(int(val["from_node"]), int(val["to_node"])) for key, val in data["valves"].items()])
    # G.add_edges_from([(int(val["from_node"]), int(val["to_node"])) for key, val in data["short_pipes"].items()])
    # G.add_edges_from([(int(val["from_node"]), int(val["to_node"])) for key, val in data["resistors"].items()])
    G.add_edges_from([(int(val["from_node"]), int(val["to_node"])) for key, val in data["control_valves"].items()])
    # G.add_edges_from([(int(val["from_node"]), int(val["to_node"])) for key, val in data["loss_resistors"].items()])
    return G

if __name__ == '__main__':
    network_data = read_data('../../data/GasLib-582/')
    G = create_graph(network_data)
    print([len(i) for i in sorted(nx.cycle_basis(G))])