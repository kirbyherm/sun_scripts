// Kirby Hermansen 
// 2022-10-19
// Purpose: To create a histogram of the sum of SuN's individual segments

#include "tas.cc"
#include "ss.cc"
#include "mult.cc"
#include "tas_ss_2d.cc"
#include "remove_daughter.cc"
#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

int main(){
    tas_parent();
    tas_daughter();
    ss_parent();
    ss_daughter();
    mult_parent();
    mult_daughter();
    tas_ss_parent();
    tas_ss_daughter();
    remove_daughter_tas();
    remove_daughter_tas_ss();
    remove_daughter_ss();
    remove_daughter_mult();
    return 0;
}
