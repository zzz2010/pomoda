/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja   *
 *   ari@ebi.ac.uk   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "intmatrix.h"

using namespace std;

IntMatrix::IntMatrix(int xa, std::string n)
{
    assert(xa>0);
    x = xa;
    y = z = w = 1;
    name = n;

    allocate();
}

IntMatrix::IntMatrix(int xa, int ya, std::string n)
{
    assert(xa>0);
    x = xa;
    assert(ya>0);
    y = ya;
    z = w = 1;
    name = n;

    allocate();
}

IntMatrix::IntMatrix(int xa, int ya, int za, std::string n)
{
    assert(xa>0);
    x = xa;
    assert(ya>0);
    y = ya;
    assert(za>0);
    z = za;
    w = 1;
    name = n;

    allocate();
}

IntMatrix::IntMatrix(int xa, int ya, int za, int wa, std::string n)
{
    assert(xa>0);
    x = xa;
    assert(ya>0);
    y = ya;
    assert(za>0);
    z = za;
    assert(wa>0);
    w = wa;
    name = n;

    allocate();
}

IntMatrix::~IntMatrix()
{
//    cout<<"int delete "<<name<<endl;
    //delete []data;
    x=y=z=w=0;
}

void IntMatrix::allocate()
{
    data = new int[x*y*z];
}

void IntMatrix::initialise(int v)
{
    FOR(i,x) {
	FOR(j,y) {
	    FOR(k,z) {
		FOR(l,w) {
		    data[i + j*x + k*x*y + l*x*y*z] = v;
		}
	    }
	}
    }
}

void IntMatrix::print(string header, string consensus, FILE *outstream)
{
    fprintf(outstream, "%s\n", header.c_str());
    //FOR(i,x) {
    for(i = 1; i < x; ++i) {
        FOR(j,y) {
            FOR(k,z) {
                FOR(l,w) {
                    fprintf(outstream, "+%i ", data[i + j*x + k*x*y + l*x*y*z]);
                }
                if(w>1)
                    fprintf(outstream, "\n");
            }
            if(z>1)
                fprintf(outstream, "\n");
        }
        if(y>1)
            fprintf(outstream, "\n");
    }
    //if(x>1) {
        //fprintf(outstream, "\n");
    //}
    
    fprintf(outstream, "#consensus:\t%s\n\n", consensus.c_str());
};


void IntMatrix::print_one_zero(string header, string consensus, FILE *outstream)
{
    fprintf(outstream, "%s\n", header.c_str());
    for(i = 1; i < x; ++i) {
        FOR(j,y) {
	    if(j == (y-1))
	    	fprintf(outstream, "|tot: ");
            FOR(k,z) {
                FOR(l,w) {
                    fprintf(outstream, "%i ", data[i + j*x + k*x*y + l*x*y*z]);
                }
                if(w>1)
                    fprintf(outstream, "\n");
            }
            if(z>1)
                fprintf(outstream, "\n");
        }
        if(y>1)
            fprintf(outstream, "\n");

    }

    fprintf(outstream, "#consensus:\t%s\n\n", consensus.c_str());
};

string IntMatrix::get_consensus(double pseudo_weight = 0, double genome_fra = 0.5, double genome_frc = 0.5 )
{
    string consensus;
    string consensus_pattern;

    double total_freq = 0;
    double matrix_value = 0;
    double matrix_weight = 0;
    double nuc_freq = 0;
    double prior_prob = 0;
    
    for(i = 1; i < x; ++i) {
	total_freq = 0;
	consensus_pattern = "";

	FOR(j,y) {
	    FOR(k,z) {
		FOR(l,w) {
		    total_freq += double(data[i + j*x + k*x*y + l*x*y*z]);
		}
	    }
	}

	FOR(j,y) {
	    nuc_freq = 0;
    	    FOR(k,z) {
		FOR(l,w) {
		    nuc_freq += double(data[i + j*x + k*x*y + l*x*y*z]);
		}
	    }

	    if ((i == 0) || (i == 3)) { prior_prob = genome_fra / 2; }
	    if ((i == 1) || (i == 2)) { prior_prob = genome_frc / 2; }
	    matrix_value = ((nuc_freq / total_freq) + (pseudo_weight * prior_prob)) / (1 + pseudo_weight);
	    matrix_weight = log(matrix_value / prior_prob);
	
	    if((matrix_weight > 0) && (j == 0)) { consensus_pattern = consensus_pattern + "A"; }
	    if((matrix_weight > 0) && (j == 1)) { consensus_pattern = consensus_pattern + "C"; }
	    if((matrix_weight > 0) && (j == 2)) { consensus_pattern = consensus_pattern + "G"; }
	    if((matrix_weight > 0) && (j == 3)) { consensus_pattern = consensus_pattern + "T"; }		    
	}
	
	if (consensus_pattern == "A" ) 		{ consensus = consensus + "A"; }
	else if (consensus_pattern == "C" )	{ consensus = consensus + "C"; }
	else if (consensus_pattern == "G" )	{ consensus = consensus + "G"; }
	else if (consensus_pattern == "T" )	{ consensus = consensus + "T"; }
	else if (consensus_pattern == "AC" )	{ consensus = consensus + "M"; }
	else if (consensus_pattern == "AG" )	{ consensus = consensus + "R"; }
	else if (consensus_pattern == "AT" )	{ consensus = consensus + "W"; }
	else if (consensus_pattern == "CG" )	{ consensus = consensus + "S"; }
	else if (consensus_pattern == "CT" )	{ consensus = consensus + "Y"; }
	else if (consensus_pattern == "GT" )	{ consensus = consensus + "K"; }
	else if (consensus_pattern == "ACG" )	{ consensus = consensus + "V"; }
	else if (consensus_pattern == "ACT" )	{ consensus = consensus + "H"; }
	else if (consensus_pattern == "AGT" )	{ consensus = consensus + "D"; }
	else if (consensus_pattern == "CGT" )	{ consensus = consensus + "B"; }
	else					{ consensus = consensus + "N"; }
    }

    return(consensus);
};

void IntMatrix::print()
{

   FOR(i,x) {
	FOR(j,y) {
	    FOR(k,z) {
		FOR(l,w) {
		    cout<<data[i + j*x + k*x*y + l*x*y*z]<<" ";
		}
		if(w>1)
		    cout<<endl;
	    }
	    if(z>1)
		cout<<endl;
	}
	if(y>1)
	    cout<<endl;
    }
    if(x>1)
	cout<<endl;
};

void IntMatrix::set_element_name(std::string this_name)
{
    element_name = this_name;
};

void IntMatrix::get_element_name(std::string& this_name)
{
    this_name = element_name;		
};
	


