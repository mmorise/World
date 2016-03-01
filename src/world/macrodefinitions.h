//-----------------------------------------------------------------------------
// Copyright (c) 2008-2011 The Department of Arts and Culture,
// The Government of the Republic of South Africa.
//
// Contributors:  Meraka Institute, CSIR, South Africa.
//                Giulio Paci.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// MODIFIED: Giulio Paci
// DATE    : 22 January 2016
//
// Some usefull macros and defines.
//-----------------------------------------------------------------------------

#ifndef WORLD_MACRODEFINITIONS_H_
#define WORLD_MACRODEFINITIONS_H_

//
// @file macrodefinitions.h
// Definitions of macros that are common to the whole World engine.
//

//
// @ingroup WorldDefines
// @defgroup WorldMacros System Macros
// Definitions of macros that are common to the whole World engine.
// @{
//

//
// Defines
//

//
// @def WORLD_BEGIN_C_DECLS
// Start block for enclosing C code for inclusion in C++ programs.
// This allows C++ programs to include the C header files of the World
// engine. @code extern "C" { @endcode
// @hideinitializer
//

//
// @def WORLD_END_C_DECLS
// End block for enclosing C code for inclusion in C++ programs.
// This allows C++ programs to include the C header files of the World
// engine. @code } @endcode
// @hideinitializer
//

#undef  WORLD_BEGIN_C_DECLS
#undef  WORLD_END_C_DECLS
#ifdef __cplusplus
#  define WORLD_BEGIN_C_DECLS extern "C" {
#  define WORLD_END_C_DECLS }
#else  // !__cplusplus
#  define WORLD_BEGIN_C_DECLS
#  define WORLD_END_C_DECLS
#endif  // __cplusplus

//
// @def WORLD_API
// @hideinitializer
// Declares a symbol to be exported for shared library usage.
//

//
// @def WORLD_LOCAL
// @hideinitializer
// Declares a symbol hidden from other libraries.
//

//
// @def WORLD_PLUGIN_API
// @hideinitializer
// Declares a symbol to be exported for plug-in usage.
//

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
#  define WORLD_HELPER_DLL_IMPORT __declspec(dllimport)
#  define WORLD_HELPER_DLL_EXPORT __declspec(dllexport)
#  define WORLD_HELPER_DLL_LOCAL
#else  // ! defined _WIN32 || defined __CYGWIN__ || defined WORLD_WIN32
#  if __GNUC__ >= 4
#    define WORLD_HELPER_DLL_IMPORT __attribute__ ((visibility("default")))
#    define WORLD_HELPER_DLL_EXPORT __attribute__ ((visibility("default")))
#    define WORLD_HELPER_DLL_LOCAL  __attribute__ ((visibility("hidden")))
#  else  // ! __GNUC__ >= 4
#    define WORLD_HELPER_DLL_IMPORT
#    define WORLD_HELPER_DLL_EXPORT
#    define WORLD_HELPER_DLL_LOCAL
#  endif  // __GNUC__ >= 4
#endif  // defined _WIN32 || defined __CYGWIN__ || defined WORLD_WIN32

//
// WORLD_LIBRARIES_EXPORTS
// ----------------------
// WORLD_LIBRARIES_EXPORTS should be defined when compiling a
// shared/dynamic library.
//
// Now we use the generic helper definitions above to define
// WORLD_API and WORLD_LOCAL. WORLD_API is used for the public API symbols.
// It's either DLL imports or DLL exports (or does nothing for static build)
// WORLD_LOCAL is used for non-api symbols.
//
// WORLD_SRC
// ------------------
// WORLD_SRC should be defined when building World (instead of just using it).
//

#ifdef WORLD_LIBRARIES_EXPORTS
#  ifdef WORLD_SRC
#    define WORLD_API WORLD_HELPER_DLL_EXPORT
#  else  // !WORLD_SRC
#    define WORLD_API WORLD_HELPER_DLL_IMPORT
#  endif  // WORLD_SRC
#  define WORLD_LOCAL WORLD_HELPER_DLL_LOCAL
#else  // !WORLD_LIBRARIES_EXPORTS (static library)
#  define WORLD_API
#  define WORLD_LOCAL
#endif  // WORLD_LIBRARIES_EXPORTS

//
// @}
// end documentation
//

#endif  // WORLD_MACRODEFINITIONS_H_
