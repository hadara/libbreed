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

function _sorted_family_tree(&$relationships) {
    # we have to ensure that no animal comes in the relationship list before its parents

    # build hashmap from animal_id to record to make sorting simpler
    $relmap = array();
    foreach ($relationships as $rel) {
        $relmap[$rel['id']] = $rel;
    }

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
        $relationships_sorted = _sorted_family_tree($relationships);
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
}
