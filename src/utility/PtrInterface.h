
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
// Framework Exception types
// Copyright(c) 1993-2006_2007, David R. Cheriton, all rights reserved.
//
// Edited by Mark Linton for CS 249A Fall 2014.
//

#ifndef PTRINTERFACE_H
#define PTRINTERFACE_H

class PtrInterface {
public:

    PtrInterface() :
        ref_(0)
    {
        // Nothing else to do.
    }

    unsigned long references() const {
        return ref_;
    }

    enum Attribute {
        nextAttributeNumber__ = 1
    };

    // DRC - support for templates
    void newRef() const {
        PtrInterface* const ptr = const_cast<PtrInterface*>(this);
        ptr->ref_ += 1;
    }

    void deleteRef() const {
        PtrInterface* const ptr = const_cast<PtrInterface*>(this);
        ptr->ref_ -= 1;
        if (ptr->ref_ == 0) {
            ptr->onZeroReferences();
        }
    }

protected:

    virtual ~PtrInterface() {
        // Nothing to do.
    }

    virtual void onZeroReferences() {
        delete this;
    }

private:

    long unsigned ref_;

};

#endif
