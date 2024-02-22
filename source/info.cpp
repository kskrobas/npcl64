/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * info.cpp
* Copyright (C) 2019 Kazimierz Skrobas <kskrobas@unipress.waw.pl>
 *
 * npcl is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * npcl is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include "info.h"
#include <iostream>
using namespace std;
void info()
{
            cerr<<" date: "<<__DATE__<<endl;
            cerr<<" author: Kazimierz Skrobas, kskrobas@unipress.waw.pl"<<endl;
            cerr<<" NPCLPATH: ";

            if(const char* env_p = std::getenv("NPCLPATH"))
               cerr<<env_p;        
}
