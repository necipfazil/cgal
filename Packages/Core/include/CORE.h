/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: CORE.h
 * Synopsis:
 *      The main inclusion file for the Core Library system.
 *      All "Core programs" must include this file.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_H
#define CORE_H

#include "CoreImpl.h"
#include "CoreAux.h"
#include "CoreDefs.h"

// User can still access machine types:
typedef double machine_double;
typedef long machine_long;

#ifndef Level 
#   define Level  DEFAULT_LEVEL
#endif

#if Level  == 1
#   define Real double
#   define Expr double
#elif Level  == 2
#   undef long
#   undef double
#   include "Real.h"
#   define long Real
#   define double Real
#   define Expr Real
#elif Level  == 3
#   undef long
#   undef double
#   include "Expr.h"
#   define long Expr
#   define double Expr
#   define Real Expr
#elif Level == 4
#   include "Expr.h"
#endif

// automaticall link necessary static library under visual c++
#ifdef _MSC_VER
  #ifdef _DEBUG
    #pragma comment(lib, "coreDebug.lib")
    #pragma comment(lib, "gmpDebug.lib")
  #else
    #pragma comment(lib, "core.lib")
    #pragma comment(lib, "gmp.lib")
  #endif
#endif

#ifndef CORE_NO_AUTOMATIC_NAMESPACE
using namespace CORE;
#endif

#endif

