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
#include "flmatrix.h"

using namespace std;

extern double sumLogs(double a, double b);

FlMatrix::FlMatrix(int xa, std::string n)
{
    assert(xa>0);
    x = xa;
    y = z = w = 1;
    name = n;

    allocate();
}

FlMatrix::FlMatrix(int xa, int ya, std::string n)
{
    assert(xa>0);
    x = xa;
    assert(ya>0);
    y = ya;
    z = w = 1;
    name = n;

    allocate();
}

FlMatrix::FlMatrix(int xa, int ya, int za, std::string n)
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

FlMatrix::FlMatrix(int xa, int ya, int za, int wa, std::string n)
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

FlMatrix::~FlMatrix()
{
//    cout<<"fl delete "<<name<<endl;
    //delete []data;
    x=y=z=w=0;
}

void FlMatrix::allocate()
{
    data = new double[x*y*z*w];
}

void FlMatrix::initialise(double v)
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

void FlMatrix::print(string header, FILE* outstream)
{
    fprintf(outstream, "%s\n", header.c_str());
    FOR(i,x) {
        FOR(j,y) {
            FOR(k,z) {
                FOR(l,w) {
                    fprintf(outstream, "+%.2f ", data[i + j*x + k*x*y + l*x*y*z]);
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
    if(x>1)
        fprintf(outstream, "\n");
};

void FlMatrix::print(string header, string consensus, FILE *outstream)
{
    fprintf(outstream, "%s\n", header.c_str());
    for(i = 1; i < x; ++i) {
        FOR(j,y) {
            FOR(k,z) {
                FOR(l,w) {
                    fprintf(outstream, "+%.2f ", data[i + j*x + k*x*y + l*x*y*z]);
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

void FlMatrix::print()
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
}

void FlMatrix::printStr(string* str)
{
    char s[10];
    FOR(i,x) {
	FOR(j,y) {
	    FOR(k,z) {
		FOR(l,w) {
		    sprintf(s,"%.1f ",data[i + j*x + k*x*y + l*x*y*z]);
		    *str += s;
		}
		if(w>1)
		    *str += "\n";
	    }
	    if(z>1)
		*str += "\n";
	}
	if(y>1)
	    *str += "\n";
    }
    if(x>1)
	*str += "\n";
}

void FlMatrix::printRevStr(string* str)
{
    string tmp;
    char s[10];
    FOR(i,x) {
	FOR(j,y) {
	    FOR(k,z) {
		FOR(l,w) {
		    sprintf(s,"%.1f ",data[i + j*x + k*x*y + l*x*y*z]);
		    tmp += s;
		}
		if(w>1)
		    tmp += "\n";
	    }
	    if(z>1)
		tmp += "\n";
	}
	if(y>1)
	    tmp += "\n";
    }
    if(x>1)
	tmp += "\n";
    *str = tmp+*str;
}

double FlMatrix::sumLogs(double a, double b) 
{

    if(a==-HUGE_VAL && b==-HUGE_VAL) {
	return -HUGE_VAL;
    }
    else if(a==-HUGE_VAL){
	return b;
    }
    else if(b==-HUGE_VAL){
	return a;
    }
    if(b>a){
	double c = a;
	a = b;
	b = c;
    }

    return (a+log(1+exp(b-a)));
}
