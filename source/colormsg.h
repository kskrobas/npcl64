/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * colormsg.h
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


#ifndef COLORMSG_H
#define COLORMSG_H

#define cred    "\u001b[48;5;1m"
#define cgreen  "\u001b[38;5;2m"
#define cyellow "\u001b[38;5;3m"
#define cblue   "\u001b[38;5;4m"
#define cdef    "\033[0m"

#include <iostream>
#include <string>

using namespace std;


void errMsg(const string msg);
void warnMsg(const string msg);
void infoMsg(const string msg);
void logMsg(const  string &msg);

#endif // TERMINALCOLORS_H
