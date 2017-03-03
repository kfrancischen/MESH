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
