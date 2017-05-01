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

#ifndef UTILITY_H
#define UTILITY_H

//
// Common framework pragmas, types, and includes.
//

#ifdef _MSC_VER

    //
    // Ignore overzealous VC++ warnings.
    //
    // 4100: Unreferenced formal parameter message.
    // 4250: Inherits symbol via dominance.
    // 4355: 'this' : used in base member initializer list
    // 4511: Copy constructor could not be generated.
    // 4512: Assignment operator could not be generated.
    // 4514: Unreferenced inline function has been removed.
    // 4522: Multiple assignment operators of a single type.
    // 4624: Destructor could not be generated because
    //       a base class destructor is inaccessible.
    // 4625: Copy constructor could not be generated because
    //       a base class copy constructor is inaccessible.
    // 4626: Assignment operator could not be generated because
    //       a base class copy constructor is inaccessible.
    //
#   pragma warning( disable : 4100 4250 4355 4511 4512 4514 4522 4624 4625 4626 )

#   define _noinline __declspec(noinline)

#   if _MSC_VER >= 1600
#       define null (nullptr)
#   else
#       define null (0)

        typedef int nullptr_t;
#   endif

#endif

#ifdef __clang__
#   define _noinline __attribute__((noinline))
#   define null (nullptr)
#endif

#ifdef __GNUC__
#   define _noinline __attribute__((noinline))
#   define null (nullptr)
#endif

#ifndef _noinline
#   define _noinline /**/
#endif

#ifndef null
#   define null (0)
#endif
#include <string.h>
#include <string>
#include <iostream>
#include <algorithm>
using std::string;
namespace UTILITY {

#   include "utility/Ptr.h"
#   include "utility/PtrInterface.h"
#   include "utility/NamedInterface.h"
#   include "utility/Exception.h"
#   include "utility/StringUtil.h"
#   include "utility/Sort.h"
}

#endif /* FWK_FWK_H */
