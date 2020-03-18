<?php

function get_item(&$array, $key, $default=null) {
    if (isset($array[$key])) {
        return $array[$key];
    }
    return $default;
}

function print_numerator_table(&$numerator_table) {
    foreach ($numerator_table as $y) {
        foreach ($y as $x) {
            echo($x."\t");
        }
        echo("\n");
    }
}

function _sorted_family_tree(&$relationships, &$relmap) {
    # we have to ensure that no animal comes in the relationship list before its parents

    # helper "set" to keep track of which animals we have already sent to the output
    $animals_seen = array();

    $relationships_out = array();
    # general idea of the sorting algorithm is to move all the animals which do not have their parents
    # present in the set to the end of the output list and delete them from the relmap
    # Do this until either relmap is empty or there hasn't been a change in relmap size during the cycle
    # in which case the remaining relationship tree contains a cycle (which means the input data is invalid)
    while (1) {
        $relmap_size_before = sizeof($relmap);
        foreach ($relmap as $animal_id => $animal_rec) {
            $parent1 = get_item($animal_rec, 'parent1');
            $parent2 = get_item($animal_rec, 'parent2');
            if (($parent1 === null || !isset($relmap[$parent1])) && ($parent2 === null || !isset($relmap[$parent2]))) {
                # if parents are referenced but not really present in the dataset then
                # we create phantom leaf nodes for them which requires less special casing in the
                # subsequent steps
                if ($parent1 !== null && !isset($animals_seen[$parent1])) {
                    array_push($relationships_out, array('id' => $parent1));
                    $animals_seen[$parent1] = true;
                }
                if ($parent2 !== null && !isset($animals_seen[$parent2])) {
                    array_push($relationships_out, array('id' => $parent2));
                    $animals_seen[$parent2] = true;
                }

                array_push($relationships_out, $animal_rec);
                $animals_seen[$animal_id] = true;
                unset($relmap[$animal_id]);
            }
        }
        if (sizeof($relmap) == 0) {
            break;
        }
        if (sizeof($relmap) == $relmap_size_before) {
            throw new Exception('relationship tree contains cycles which should not be possible!');
        }
    }
    return $relationships_out;
}

class Pedigree {
    function __construct(&$relationships) {
        # build hashmap from animal_id to record to make sorting simpler
        $relmap = array();
        foreach ($relationships as $rel) {
            $relmap[$rel['id']] = $rel;
        }
        $this->relmap = $relmap;

        $relationships_sorted = _sorted_family_tree($relationships, $relmap);
        $this->numerator_table = $this->get_numerator_table($relationships_sorted);
    }

    function parent_numerator($animal1, $animal2, &$numerator_table) {
        if (isset($numerator_table[$animal1][$animal2])) {
            return $numerator_table[$animal1][$animal2];
        }
        return 0;
    }

    function get_numerator_table(&$relationships_sorted) {
        // initialize table, not strictly necessary but makes subsequent code easier
        $numerator_table = array();
        foreach ($relationships_sorted as $rel_y) {
            $numerator_table[$rel_y['id']] = array();
            foreach ($relationships_sorted as $rel_x) {
                $numerator_table[$rel_y['id']][$rel_x['id']] = null;
            }
        }
    
        foreach ($relationships_sorted as $rel_y) {
            foreach ($relationships_sorted as $rel_x) {
                if ($numerator_table[$rel_y['id']][$rel_x['id']] !== null) {
                    // entries under the diagonal are already filles
                    continue;
                }
                if ($rel_x['id'] === $rel_y['id']) {
                    $numerator = 1 + 0.5*$this->parent_numerator(get_item($rel_x, 'parent2'), get_item($rel_x, 'parent1'), $numerator_table);
                } else {
                    // lets check how y is related to parents on x
                    $parent1_numerator = 0.5 * $this->parent_numerator($rel_y['id'], get_item($rel_x, 'parent2'), $numerator_table);
                    $parent2_numerator = 0.5* $this->parent_numerator($rel_y['id'], get_item($rel_x, 'parent1'), $numerator_table);
                    $numerator = $parent1_numerator + $parent2_numerator;
                }
                $numerator_table[$rel_y['id']][$rel_x['id']] = $numerator;
                $numerator_table[$rel_x['id']][$rel_y['id']] = $numerator;
            }
        }
        return $numerator_table;
    }

    function array_intersect_closest(&$array1, &$array2) {
        // like std. array_intersect_key() but keeps records with lowest value
        $retval = array();
        foreach ($array1 as $k1 => $v1) {
            if (!isset($array2[$k1])) {
                continue;
            }

            if ($v1 < $array2[$k1]) {
                $retval[$k1] = $array1[$k1];
            } else {
                $retval[$k1] = $v1;
            }

        }
        return $retval;
    }

    function find_all_ancestors($animal_id, $depth=0) {
        $animal_rec = get_item($this->relmap, $animal_id);
        if ($animal_rec === null) {
            // we might be called on a dataset where phantom leaf nodes aren't added for parents that
            // do not exist in the dataset
            return array($animal_id => $depth);
        }

        if (isset($animal_rec['parent1'])) {
            $parent1_ancestors = $this->find_all_ancestors($animal_rec['parent1'], $depth+1);
        } else {
            $parent1_ancestors = array();
        }

        if (isset($animal_rec['parent2'])) {
            $parent2_ancestors = $this->find_all_ancestors($animal_rec['parent2'], $depth+1);
        } else {
            $parent2_ancestors = array();
        }

        // in the case where same parent is seen on both sides of the family tree we only keep the one that is closest
        $shared_ancestors = $this->array_intersect_closest($parent1_ancestors, $parent2_ancestors);
        foreach ($shared_ancestors as $k => $v) {
            if ($parent1_ancestors[$k] <= $parent2_ancestors[$k]) {
                unset($parent2_ancestors[$k]);
            } else {
                unset($parent1_ancestors[$k]);
            }
        }
        if ($depth > 0) {
            $self = array($animal_id => $depth);
        } else {
            $self = array();
        }
        $retval = $self + $parent1_ancestors + $parent2_ancestors;
        return $retval;
    }

    // interface
    function inbreeding_coefficient($animal_id) {
        return $this->numerator_table[$animal_id][$animal_id] - 1;
    }
    
    function wright_coefficient_of_relationship($animal1_id, $animal2_id) {
        return $this->numerator_table[$animal1_id][$animal2_id] / 
                (sqrt(1+$this->inbreeding_coefficient($animal1_id)) * 
                 sqrt(1+$this->inbreeding_coefficient($animal2_id))
                );
    }

    function common_ancestors($animal1_id, $animal2_id) {
        /* returns common ancestors for given animals as a list of common ancestors ordered
         * from closest to farthest. Each common ancestor is represented by a dict in the form of:
         *  {'id' => <animal_id>, 'distance': <integer indicating minimal edge count in the graph to closest common ancestor>}
         * Empty list is returned if no common ancestors are found
         */
        $animal1_ancestors = $this->find_all_ancestors($animal1_id);
        $animal2_ancestors = $this->find_all_ancestors($animal2_id);
        $shared_ancestors = $this->array_intersect_closest($animal1_ancestors, $animal2_ancestors);
        asort($shared_ancestors, SORT_NUMERIC);
        $retval = array();
        foreach ($shared_ancestors as $key => $val) {
            array_push($retval, array('id' => $key, 'distance' => $val));
        }
        return $retval;
    }

    function most_recent_common_ancestor($animal1_id, $animal2_id) {
        /* returns ID of the most recent common ancestor for given animals
         * null if unknown.
         * When multiple ancestors have the same distance then it's undefined which one of these is returned
         */
        $shared_ancestors = $this->common_ancestors($animal1_id, $animal2_id);
        if (count($shared_ancestors) > 0) {
            return $shared_ancestors[0]['id'];
        }
        return null;
    }
}
