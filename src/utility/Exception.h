
/* This file is part of MESH
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

#ifndef EXCEPTION_H
#define EXCEPTION_H

class Exception {
public:

    string what() const {
        return what_;
    }

    //
    // Providing a virtual destructor enables us to catch an Exception, and
    // then dynamic_cast it to a derived exception type or fetch its typeid, etc.
    //
    virtual ~Exception() {
        // Nothing else to do.
    }

protected:

    Exception(char const * str) :
        what_(str)
    {
        // Nothing else to do.
    }

    Exception(const string& str) :
        what_(str)
    {
        // Nothing else to do.
    }

    Exception(const Exception& e) :
        what_(e.what_)
    {
        // Nothing else to do.
    }

private:

    string what_;

};


std::ostream& operator<<(std::ostream&, Exception const &);


class UnknownTypeException : public Exception {
public:

    UnknownTypeException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class UnknownAttrException : public Exception {
public:

    UnknownAttrException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class UnknownDelimiterException : public Exception {
public:

    UnknownDelimiterException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class UnknownArgException : public Exception {
public:

    UnknownArgException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class InternalException : public Exception {
public:

    InternalException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class RangeException : public Exception {
public:

    RangeException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class MemoryException : public Exception {
public:

    MemoryException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class StorageException : public Exception {
public:

    StorageException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class TimeoutException : public Exception {
public:

    TimeoutException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class NameInUseException : public Exception {
public:

    NameInUseException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class IllegalNameException : public Exception {
public:

    IllegalNameException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class PermissionException : public Exception {
public:

    PermissionException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class NoImplementationException : public Exception {
public:

    NoImplementationException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};


class AttributeNotSupportedException : public NoImplementationException {
public:

    // python attempts to access unimplemented attributes.
    AttributeNotSupportedException(const string& what) :
        NoImplementationException(what)
    {
        // Nothing else to do.
    }

};


class NullPointerException : public Exception {
public:

    NullPointerException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};

class FileNotExistException : public Exception {
public:

    FileNotExistException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};

class ValueException : public Exception {
public:

    ValueException(const string& what) :
        Exception(what)
    {
        // Nothing else to do.
    }

};
static void checkNull(void* const ptr) {
    if (ptr == null) {
        throw NullPointerException("Null pointer");
    }
}

static void checkNull(const void* const ptr) {
    if (ptr == null) {
        throw NullPointerException("Null pointer");
    }
}


#endif /* EXCEPTION_H */
