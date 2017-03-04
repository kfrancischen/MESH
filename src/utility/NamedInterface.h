
/* Copyright (C) 2016-2018, Stanford University
 * This file is part of MESH
 *
 * MESH is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * MESH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
 //  this code is adopted from CS 249A
// Copyright (c) 1993-2007 David R. Cheriton, all rights reserved.
//
// Edited by Mark Linton for CS 249A Fall 2014.
//

#ifndef FWK_NAMEDINTERFACE_H
#define FWK_NAMEDINTERFACE_H

/**
 * NamedInterface is an abstract base class for objects that have a string name
 * initialized during object construction.
 */
class NamedInterface : public PtrInterface {
public:

    const string& name() const {
        return name_;
    }

protected:

    NamedInterface(const string& name) :
        name_(name)
    {
        // Nothing else to do.
    }

private:

    string name_;

};

#endif
