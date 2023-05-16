/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cprogess.h
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


#ifndef CPROGESS_H
#define CPROGESS_H


#include <thread>
#include <iostream>
#include <functional>
#include <string>



class CProgress final
{
private:
        std::thread  pthread;
        volatile static float pos;
        volatile static bool active;

        float istep=1;
        
        static bool printShown;
        static void print();

public:
    //cprogess();
    //ThreadWrapper(const ThreadWrapper&) = delete;


        CProgress();
        CProgress(const CProgress & ) = delete;
        CProgress & operator=(const CProgress &) = delete;

        CProgress(CProgress &&obj);
        CProgress & operator=(CProgress && );

        ~CProgress();
        //CProgress(std::function<void()> func);



        void start(const size_t numOfsteps);
        void stop();
        void next();


        void   operator++ ();
        void   operator++ (int );

        static std::string linehead;
        static std::string title;
        static size_t stripeWidth;


};

#endif // CPROGESS_H
