# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 11:59:22 2015
@author: Ingo Scholtes

(c) Copyright ETH Zurich, Chair of Systems Design, 2015-2017
"""
import pathpy as pp
import pytest


slow = pytest.mark.slow


def test_readfile_import(path_from_ngram_file):
    levels = list(path_from_ngram_file.paths.keys())
    max_level = max(levels)
    expected = 5
    assert max_level == expected, \
        "The nodes have not been imported correctly"

    assert path_from_ngram_file.getNodes() == {'a', 'b', 'c', 'd', 'e'}, \
        "Wrong node labels"


def test_write_file(tmpdir, random_paths):
    dir_path = tmpdir.mkdir("sub").join("test.edges")
    p = random_paths(30, 50)

    expected_seq = ''.join(p.getSequence())
    expected_paths = sorted(expected_seq.split('|'))

    p.writeFile(dir_path.strpath)
    p2 = pp.Paths.readFile(dir_path.strpath, pathFrequency=True)

    read_back = ''.join(p2.getSequence())
    read_back_paths = sorted(read_back.split('|'))

    assert expected_paths == read_back_paths


def test_read_edges_import(path_from_edge_file):
    """test if the Paths.readEdges functions works"""
    levels = list(path_from_edge_file.paths.keys())
    max_level = max(levels)
    expected_max = 1
    assert expected_max == max_level, \
        "The nodes have not been imported correctly"

    assert path_from_edge_file.getNodes() == {'1', '2', '3', '5'}, \
        "Nodes not imported correctly"


def test_read_edges_undirected(path_from_edge_file_undirected):
    p = path_from_edge_file_undirected
    layers = list(p.paths.keys())
    max_layers = max(layers)
    expected_layers = 1
    assert max_layers == expected_layers, \
        "The nodes have not been imported correctly"

    assert p.getNodes() == {'1', '2', '3', '5'}, \
        "Nodes not imported correctly"


def test_get_sequence(path_from_ngram_file):
    from collections import Counter
    """Test if the Paths.getSequence function works correctly"""
    sequence = path_from_ngram_file.getSequence()
    sequence = "".join(sequence)
    ct = Counter(sequence.split('|'))
    assert dict(ct) == {'': 1, 'abcdab': 2, 'dedab': 4}, \
        "Returned the wrong sequence"


def test_get_unique_paths(random_paths):
    p = random_paths(90, 90)
    assert p.getUniquePaths() == 87, \
        "Wrong number of paths detected"


def test_observation_count_file(path_from_ngram_file):
    assert path_from_ngram_file.ObservationCount() == 6, \
        "Wrong number of observations detected"


def test_observation_count_large(random_paths):
    p = random_paths(90, 90)
    assert p.ObservationCount() == 193, \
        "Wrong number of observations detected"


def test_path_summary(random_paths):
    p = random_paths(90, 90)
    print(p)


def test_summary_multi_order_model(random_paths):
    p = random_paths(90, 90)
    multi = pp.MultiOrderModel(paths=p, maxOrder=3)
    print(multi)


def test_get_shortest_paths(path_from_ngram_file):
    path_from_ngram_file.getShortestPaths()
    paths_dict = path_from_ngram_file.getShortestPaths()
    expected_paths = {('d', 'a'): {('d', 'a')},
                      ('b', 'd'): {('b', 'c', 'd')},
                      ('d', 'e'): {('d', 'e')},
                      ('a', 'c'): {('a', 'b', 'c')},
                      ('a', 'a'): {('a',)},
                      ('e', 'a'): {('e', 'd', 'a')},
                      ('e', 'b'): {('e', 'd', 'a', 'b')},
                      ('e', 'e'): {('e',)},
                      ('a', 'b'): {('a', 'b')},
                      ('b', 'b'): {('b',)},
                      ('c', 'd'): {('c', 'd')},
                      ('d', 'b'): {('d', 'a', 'b')},
                      ('c', 'a'): {('c', 'd', 'a')},
                      ('b', 'a'): {('b', 'c', 'd', 'a')},
                      ('c', 'b'): {('c', 'd', 'a', 'b')},
                      ('e', 'd'): {('e', 'd')},
                      ('a', 'd'): {('a', 'b', 'c', 'd')},
                      ('d', 'd'): {('d',)},
                      ('c', 'c'): {('c',)},
                      ('b', 'c'): {('b', 'c')}
                      }
    paths_to_check = dict()
    for k in paths_dict:
        for p in paths_dict[k]:
            paths_to_check[(k, p)] = paths_dict[k][p]
    assert paths_to_check == expected_paths


def test_get_contained_paths():
    path_to_check = ('a', 'b', 'c', 'd', 'e', 'f', 'g')
    node_filter = ('a', 'b', 'd', 'f', 'g')
    cont_paths = pp.Paths.getContainedPaths(path_to_check, node_filter)
    expected = [('a', 'b'), ('d',), ('f', 'g')]
    assert cont_paths == expected


def test_filter_paths(path_from_ngram_file):
    from collections import Counter
    p = path_from_ngram_file
    new_paths = p.filterPaths(node_filter=['a', 'b', 'c'])
    expected_sequence = {'': 1, 'ab': 6, 'abc': 2}

    new_sequence = ''.join(new_paths.getSequence())
    ct = Counter(new_sequence.split('|'))
    assert dict(ct) == expected_sequence


def test_project_paths(path_from_ngram_file):
    from collections import Counter
    p = path_from_ngram_file
    mapping = {'a': 'x', 'b': 'x', 'c': 'y', 'd': 'y', 'e': 'y'}
    new_p = p.projectPaths(mapping=mapping)
    new_sequence = ''.join(new_p.getSequence())
    ct = Counter(new_sequence.split('|'))
    expected_sequence = {'': 1, 'xxyyxx': 2, 'yyyxx': 4}
    assert dict(ct) == expected_sequence


def test_get_nodes(random_paths):
    p = random_paths(3, 9)
    rest = p.getNodes()
    expected = {'b', 'o', 'u', 'v', 'w', 'y'}
    assert rest == expected
