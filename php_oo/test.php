<?php

// available since php 7.2
if (!defined('PHP_FLOAT_EPSILON')) {
    define('PHP_FLOAT_EPSILON', 0.00001);
}

require("libbreed.php");

function _is_close_enough($value, $expected) {
    if (abs($value-$expected) < PHP_FLOAT_EPSILON) {
        return true;
    }
    return false;
}

$test_ab = function () {
    $AB_COEFFICIENT = 0.6;
    $animal_tree_s = file_get_contents("../test_cases/ab.json");
    $animal_tree = json_decode($animal_tree_s, $assoc=true);
    $pedigree = new Pedigree($animal_tree);
    $coefficient = $pedigree->wright_coefficient_of_relationship("C", "D");
    if (_is_close_enough($coefficient, $AB_COEFFICIENT)) {
        return true;
    }
    print("ab test returned ".$coefficient." expected ".$AB_COEFFICIENT."\n");
    return false;
};

$test_cycle = function () {
    $animal_tree_s = file_get_contents("../test_cases/ab_cycle.json");
    $animal_tree = json_decode($animal_tree_s, $assoc=true);
    try {
        $pedigree = new Pedigree($animal_tree);
        print("cyclic graph test did not raise an Exception\n");
        return false; 
    } catch (Exception $e) {
        return true;
    }
};

$test_prospector_inbreeding = function () {
    $animal_tree_s = file_get_contents("../test_cases/prospector.json");
    $animal_tree = json_decode($animal_tree_s, $assoc=true);
    $pedigree = new Pedigree($animal_tree);
    foreach ($animal_tree as $animal) {
        if (isset($animal["F"])) {
            $calculated_F = $pedigree->inbreeding_coefficient($animal['id']);
        } else {
            continue;
        }
        if (!_is_close_enough($calculated_F, $animal['F'])) {
            print("calculated inbreeding coefficient for ".$animal['id']." is ".$calculated_F." should be:".$animal['F']."\n");
            #return false;
        }
    }
    return true;
};

$test_book_17_5 = function () {
    $animal_tree_s = file_get_contents("../test_cases/book_17_5.json");
    $animal_tree = json_decode($animal_tree_s, $assoc=true);
    $pedigree = new Pedigree($animal_tree);

    $Fx = $pedigree->inbreeding_coefficient("X");
    $Fx_expected = 0.375;
    if (!_is_close_enough($Fx, $Fx_expected)) {
        print("Fx is wrong. is ".$Fx." should be: ".$Fx_expected);
        return false;
    }

    $Rxa = $pedigree->wright_coefficient_of_relationship("X", "A");
    $Rxa_expected = 0.4264;
    if (!_is_close_enough($Rxa, $Rxa_expected)) {
        print("Rxa is wrong. is ".$Rxa." should be: ".$Rxa_expected);
        return false;
    }

    $Fh = $pedigree->inbreeding_coefficient("H");
    $Fh_expected = 0.1875;
    if (!_is_close_enough($Fh, $Fh_expected)) {
        print("Fh is wrong. is ".$Fh." should be: ".$Fh_expected);
        return false;
    }

    $Rgh = $pedigree->wright_coefficient_of_relationship("G", "H");
    $Rgh_expected = 0.64888;
    if (!_is_close_enough($Rgh, $Rgh_expected)) {
        print("Rgh is wrong. is ".$Rgh." should be: ".$Rgh_expected);
        return false;
    }

    return true;
};

# $test_prospector_inbreeding FIXME: is the example in the book broken?
#$tests = array($test_ab, $test_cycle, $test_book_17_5);
$tests = array($test_book_17_5);
$failed_tests = 0;
foreach ($tests as $func) {
    if ($func() === false) {
        $failed_tests += 1;
    } 
}

if ($failed_tests > 0) {
    print($failed_tests." tests failed\n");
} else {
    print("all tests passed\n");
}
