import pathpy as pp
import pytest
import numpy as np
import os

test_directory = os.path.dirname(os.path.abspath(__file__))
test_data_dir = os.path.join(test_directory, 'test_data')


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",
        help="run slow tests")


def pytest_runtest_setup(item):
    if 'slow' in item.keywords and not item.config.getvalue("runslow"):
        pytest.skip("need --runslow option to run")



@pytest.fixture()
def test_data_directory():
    return test_data_dir


@pytest.fixture()
def path_from_ngram_file():
    """load the example file as pypath.Path"""
    ngram_file_path = os.path.join(test_data_dir, 'ngram_simple.ngram')
    path = pp.Paths.readFile(ngram_file_path, pathFrequency=True)
    return path


@pytest.fixture()
def path_from_edge_file():
    file_path = os.path.join(test_data_dir, 'edge_frequency.edge')
    path = pp.Paths.readEdges(file_path, weight=True)

    return path


@pytest.fixture()
def path_from_edge_file_undirected():
    file_path = os.path.join(test_data_dir, 'edge_frequency.edge')
    path = pp.Paths.readEdges(file_path, weight=True, undirected=True)
    return path


def generate_random_path(size, rnd_seed):
    """Generate a Path with random path sequences"""
    import string
    node_set = string.ascii_lowercase

    def random_ngram(p_len, nodes):
        num_elements = len(nodes)
        sequence = np.random.choice(num_elements, p_len)
        path = [nodes[i] for i in sequence]
        return ','.join(path)

    np.random.seed(rnd_seed)
    paths = pp.Paths()
    for _ in range(size):
        frequency = np.random.randint(1, 4)
        path_length = np.random.randint(1, 10)
        path_to_add = random_ngram(path_length, node_set)
        paths.addPath(path_to_add, pathFrequency=frequency)

    return paths


@pytest.fixture(scope='function')
def random_paths():
    """Generate a Path with random path sequences"""
    return generate_random_path


@pytest.fixture()
def temporal_network_object():
    t = pp.TemporalNetwork()
    # Path of length two
    t.addEdge("c", "e", 1)
    t.addEdge("e", "f", 2)

    # Path of length two
    t.addEdge("a", "e", 3)
    t.addEdge("e", "g", 4)

    # Path of length two
    t.addEdge("c", "e", 5)
    t.addEdge("e", "f", 6)

    # Path of length two
    t.addEdge("a", "e", 7)
    t.addEdge("e", "g", 8)

    # Path of length two
    t.addEdge("c", "e", 9)
    t.addEdge("e", "f", 10)

    # The next two edges continue the previous path to ( c-> e-> f-> e -> b )
    t.addEdge("f", "e", 11)
    t.addEdge("e", "b", 12)

    # This is an isolated edge (i.e. path of length one)
    t.addEdge("e", "b", 13)

    # Path of length two
    t.addEdge("c", "e", 14)
    t.addEdge("e", "f", 15)

    # Path of length two
    t.addEdge("b", "e", 16)
    t.addEdge("e", "g", 17)

    # Path of length two
    t.addEdge("c", "e", 18)
    t.addEdge("e", "f", 19)

    # Path of length two
    t.addEdge("c", "e", 20)
    t.addEdge("e", "f", 21)

    return t
