/* 
 * Scythe Statistical Library Copyright (C) 2000-2002 Andrew D. Martin
 * and Kevin M. Quinn; 2002-present Andrew D. Martin, Kevin M. Quinn,
 * and Daniel Pemstein.  All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify under the terms of the GNU General Public License as
 * published by Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.  See the text files
 * COPYING and LICENSE, distributed with this source code, for further
 * information.
 * --------------------------------------------------------------------
 *  scythestat/error.h
 */
 
/*! \file error.h 
 *
 * \brief Definitions of Scythe exception classes.
 *
 * This file contains the class definitions for
 * scythe::scythe_exception and its children.  These exception classes
 * describe all of the error conditions generated by Scythe library
 * routines.
 *
 * Furthermore, error.h contains a series of macro definitions that
 * regulate the inclusion of the library's error checking code in
 * compiled code.  These macros are controlled by the compiler flag
 * SCYTHE_DEBUG and define four levels of scythe debug
 * info, SCYTHE_DEBUG = 0, 1, 2, or 3.  The library uses these macros to
 * specify the debug level of thrown exceptions.  If we are at level
 * three, all throws are expanded into actual code, at level 2 only
 * SCYTHE_THROW_10 AND SCYTHE_THROW_20 calls are expanded, and so on.
 * Scythe developers should balance exception importance and
 * efficiency costs when making exception level choices.  For example,
 * bounds checking in matrices is done at level three primarily
 * because the added branch results in high performance penalties and
 * out-of-bounds errors shouldn't occur in well-written code, while
 * conformance checks in matrix multiplication are level 1 because the
 * checks result in little overhead relative to the cost of matrix
 * multiplication and conformation errors are easy to introduce by
 * accident.  At level 0, the library performs virtually no error
 * checking.
 *
 * While the various SCYTHE_THROW, SCYTHE_CHECK, and SCYTHE_WARN
 * macros will only typically be used by library developers, users
 * should make extensive use the tiered error reporting in Scythe by
 * setting the compiler flag SCYTHE_DEBUG.  If not explicitly set by
 * the user, the SCYTHE_DEBUG level is automatically set to 3.
 */

#ifndef SCYTHE_ERROR_H
#define SCYTHE_ERROR_H

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

/*! @cond */
#ifdef SCYTHE_DEBUG_LIB
#define SCYTHE_DEBUG_MSG(MSG)                             \
{ std::cout << "SCYTHE_DEBUG_LIB: " << MSG << std::endl; }
#else
#define SCYTHE_DEBUG_MSG(MSG)
#endif
/*! @endcond */

#define SCYTHE_THROW(EXCEP,MSG)                           \
	{                                                       \
		std::stringstream _SCYTHE_DEBUG_ss;                   \
		_SCYTHE_DEBUG_ss << MSG;                              \
		throw EXCEP(__FILE__, __func__, __LINE__,  \
				_SCYTHE_DEBUG_ss.str());                          \
	}

#define SCYTHE_CHECK(CHECK,EXCEP,MSG)                     \
{                                                         \
	if (CHECK)                                              \
		SCYTHE_THROW(EXCEP,MSG)                               \
}

/*! @cond */
#ifndef SCYTHE_DEBUG
#define SCYTHE_DEBUG 3
#endif
/*! @endcond */

#if SCYTHE_DEBUG > 0
#define SCYTHE_CHECK_10(CHECK,EXCEP,MSG) SCYTHE_CHECK(CHECK,EXCEP,MSG)
#else
#define SCYTHE_CHECK_10(CHECK, EXCEP, MSG)
#endif

#if SCYTHE_DEBUG > 1
#define SCYTHE_CHECK_20(CHECK,EXCEP,MSG) SCYTHE_CHECK(CHECK,EXCEP,MSG)
#else
#define SCYTHE_CHECK_20(CHECK, EXCEP, MSG)
#endif

#if SCYTHE_DEBUG > 2
#define SCYTHE_CHECK_30(CHECK,EXCEP,MSG) SCYTHE_CHECK(CHECK,EXCEP,MSG)
#else
#define SCYTHE_CHECK_30(CHECK, EXCEP, MSG)
#endif

#if SCYTHE_DEBUG > 0 
#define SCYTHE_THROW_10(EXCEP,MSG) SCYTHE_THROW(EXCEP,MSG)
#else
#define SCYTHE_THROW_10(EXCEP,MSG)
#endif

#if SCYTHE_DEBUG > 1 
#define SCYTHE_THROW_20(EXCEP,MSG) SCYTHE_THROW(EXCEP,MSG)
#else
#define SCYTHE_THROW_20(EXCEP,MSG)
#endif

#if SCYTHE_DEBUG > 2 
#define SCYTHE_THROW_30(EXCEP,MSG) SCYTHE_THROW(EXCEP,MSG)
#else
#define SCYTHE_THROW_30(EXCEP,MSG)
#endif

#define SCYTHE_WARN(MSG)                                              \
  {                                                                   \
  std::cerr << "WARNING in " << __FILE__ << ", "                      \
    << __func__ << ", " << __LINE__ << ": "                \
    << MSG << "\n";                                                   \
  }

#define SCYTHE_CHECK_WARN(CHECK,MSG)                                  \
  {                                                                   \
  if (CHECK)                                                          \
    SCYTHE_WARN(MSG)                                                  \
  }


namespace scythe
{
	/* Forward declaration for serr */
	class scythe_exception;

  /**** This file-local variable holds the output of the last
   * scythe_exception constructed.
   ****/
#ifdef __MINGW32__
  static scythe_exception *serr;
#else
  namespace
  {
    scythe_exception *serr;
  }
#endif

  /**** A replacement for the default terminate handler.  This outputs
   * the string held in serr before calling abort, thereby notifying
   * the user of why the program crashed.
   ****/
  inline void scythe_terminate ();

  /**** The scythe exception abstract base class ****/
  /*!
   * \brief The Scythe exception abstract base class.
   *
   * The is the base class in Scythe's error handling class tree.
   * This class extends std::exception and provides fields for
   * information about the exception, including where the exception
   * occurred in the library and a message describing the error.
   */
  class scythe_exception:public std::exception
  {
  public:
    scythe_exception (const std::string & head,
                      const std::string & file,
                      const std::string & function,
                      const unsigned int &line,
                      const std::string & message = "",
                      const bool & halt = false) throw ()
      : exception (),
        head_ (head),
        file_ (file),
        function_ (function),
        line_ (line),
        message_ (message),
				call_files_ (),
				call_funcs_ (),
				call_lines_ ()
    {
      std::ostringstream os;
      os << head_ << " in " << file_ << ", " << function_ << ", "
        << line_ << ": " << message_ << "!\n\n";

			serr = this;
      std::set_terminate (scythe_terminate);
      if (halt)
        std::terminate ();
    }

    scythe_exception (const scythe_exception & e) throw ()
      : exception (),
        head_ (e.head_),
        file_ (e.file_),
        function_ (e.function_),
        line_ (e.line_),
        message_ (e.message_),
				call_files_ (e.call_files_),
				call_funcs_ (e.call_funcs_),
				call_lines_ (e.call_lines_)
    {
    }

    scythe_exception & operator= (const scythe_exception & e) throw ()
    {
      head_ = e.head_;
      file_ = e.file_;
      function_ = e.function_;
      line_ = e.line_;
      message_ = e.message_;

      return *this;
    }

    virtual ~ scythe_exception () throw ()
    {
    }

    virtual const char *what () const throw ()
    {
      std::ostringstream os;
			for (int i = call_files_.size() - 1; i > -1; ++i) {
				os << "Called from " << call_files_[i] << ", "
					<< call_funcs_[i] << ", " << call_lines_[i] << std::endl;
			}
      os << head_ << " in " << file_ << ", " << function_ << ", "
        << line_ << ": " << message_ << "!";
      return os.str ().c_str ();
    }

    virtual std::string message () const throw ()
    {
      return message_;
    }

		virtual void add_caller (const std::string &file,
			const std::string &function, const unsigned int &line) throw ()
		{

			/* This if allows one to catch and rethrow an error in the same
			 * function w/out messing things up.  Nice to keep try-catch
			 * blocks to a minimum
			 */

			if (file != file_ && function != function_) {
				call_files_.push_back(file);
				call_funcs_.push_back(function);
				call_lines_.push_back(line);
			}
		}

  private:
    std::string head_;
    std::string file_;
    std::string function_;
    unsigned int line_;
    std::string message_;
		std::vector<std::string> call_files_;
		std::vector<std::string> call_funcs_;
		std::vector<unsigned int> call_lines_;
  };


  /**** Exception class types, added as needed ****/

  /*! 
   * \brief Memory allocation error.
   *
   *  Library members throw this exception in response to insufficient
   *  memory conditions, such as when one attempts to create a Matrix
   *  object that is bigger than available memory.
   */
  class scythe_alloc_error:public scythe_exception
  {
  public:
    scythe_alloc_error (const std::string & file,
                        const std::string & function,
                        const unsigned int &line,
                        const std::string & message = "",
                        const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE_ALLOCATION_ERROR", file, function,
          line, message, halt)
    {
    }
  };

  /*! 
   * \brief Invalid function argument.
   *
   * Library members throw this exception when callers pass incorrect
   * arguments to a function, such as when one calls the factorial
   * method with an argument less than 0.
   */
  class scythe_invalid_arg:public scythe_exception
  {
  public:
    scythe_invalid_arg (const std::string & file,
                        const std::string & function,
                        const unsigned int &line,
                        const std::string & message = "",
                        const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE_INVALID ARGUMENT", file, function,
          line, message, halt)
    {
    }
  };

  /*! 
   * \brief File i/o error.
   *
   * Library members throw this exception when errors occur during
   * file reading, writing, or creation, such as when one passes an
   * invalid file name to the Matrix class's save method.
   */
  class scythe_file_error:public scythe_exception
  {
  public:
    scythe_file_error(const std::string & file,
                       const std::string & function,
                       const unsigned int &line,
                       const std::string & message = "",
                       const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE FILE ERROR", file, function, line, 
          message, halt)
    {
    }
  };

  /*! \brief Matrix conformation error.
   *
   * Library members throw this exception when a caller passes
   * non-conforming matrices (matrices of incompatible dimensions) to
   * a routine, such as when one attempt two row vectors.
   */
  class scythe_conformation_error:public scythe_exception
  {
  public:
    scythe_conformation_error(const std::string & file,
                               const std::string & function,
                               const unsigned int &line,
                               const std::string & message = "",
                               const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE CONFORMATION ERROR", file, function, 
          line, message, halt)
    {
    }
  };

  /*! \brief Matrix dimension error.
   *
   * Library members throw this exception when a caller passes a
   * Matrix of the wrong size or shape to a routine.  For example,
   * trying to take the Cholesky decomposition of a non-square Matrix
   * causes this error.
   */

  class scythe_dimension_error:public scythe_exception
  {
  public:
    scythe_dimension_error (const std::string & file,
                            const std::string & function,
                            const unsigned int &line,
                            const std::string & message = "",
                            const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE DIMENSION ERROR", file, function,
          line, message, halt)
    {
    }
  };

  /*! \brief Null Matrix error.
   *
   * Library members throw this exception when a caller passes a null
   * Matrix to a routine when it expects a non-null argument.  For
   * example, taking the inverse of a null Matrix is impossible,
   * resulting in this exception.
   */
  class scythe_null_error:public scythe_exception
  {
  public:
    scythe_null_error(const std::string & file,
                       const std::string & function,
                       const unsigned int &line,
                       const std::string & message = "",
                       const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE NULL ERROR", file, function, line,
          message, halt)
    {
    }
  };

  /*! \brief Matrix type error.
   *
   * Library members throw this exception when a caller passes a
   * Matrix that does not satisfy some required property to a routine.
   * For example, Cholesky decomposition is designed to work on
   * positive definite matrices; trying to perform Cholesky
   * decomposition on a Matrix that does not satisfy this requirement
   * causes this exception.
   */
  class scythe_type_error:public scythe_exception
  {
  public:
    scythe_type_error(const std::string & file,
                       const std::string & function,
                       const unsigned int &line,
                       const std::string & message = "",
                       const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE TYPE ERROR", file, function, line,
          message, halt)
    {
    }
  };

  /*! \brief Element out of bounds error.
   *
   * Library members throw this exception when a caller attempts to
   * access an element outside the bounds of a data structure, such as
   * when one tries to access the 1000th element of a 200-element
   * Matrix.
   */
  class scythe_bounds_error:public scythe_exception
  {
  public:
    scythe_bounds_error(const std::string & file,
                               const std::string & function,
                               const unsigned int &line,
                               const std::string & message = "",
                               const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE BOUNDS ERROR", file, function,
          line, message, halt)
    {
    }
  };

  /*! \brief Numerical convergence error.
   *
   * Library members throw this exception when a numerical algorithm
   * fails to converge to a stable value.  For example, the BFGS
   * optimization routine throws this exception when it cannot locate
   * the minimum of a function to a given tolerance.
   */
  class scythe_convergence_error:public scythe_exception
  {
  public:
    scythe_convergence_error (const std::string & file,
                              const std::string & function,
                              const unsigned int &line,
                              const std::string & message = "",
                              const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE CONVERGENCE ERROR", file, function,
          line, message, halt)
    {
    }
  };

  /*! \brief Numerical underflow or overflow error.
   *
   * Library members throw this exception when the result of a
   * calculation, assignment, or other operation is to small or large
   * for the data type holding the value.  For example, passing
   * certain values to the gammafn function can result in underflow or
   * overflow conditions in the resulting calculations.
   */
  class scythe_range_error:public scythe_exception
  {
  public:
    scythe_range_error (const std::string & file,
                        const std::string & function,
                        const unsigned int &line,
                        const std::string & message = "",
                        const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE RANGE ERROR", file, function, line,
          message, halt)
    {
    }
  };

  /*! \brief Numerical precision error.
   *
   * Library members throw this exception when a routine cannot
   * complete a computation effectively and will sacrifice reasonable
   * precision as a consequence.  For example, passing a value too
   * close to a negative integer to the gammafn function renders the
   * function incapable of returning an accurate result and thus
   * generates this exception.
   */
  class scythe_precision_error:public scythe_exception
  {
  public:
    scythe_precision_error (const std::string & file,
                            const std::string & function,
                            const unsigned int &line,
                            const std::string & message = "",
                            const bool & halt = false) throw ()
      : scythe_exception ("SCYTHE PRECISION ERROR", file, function,
          line, message, halt)
    {
    }
  };
  
  /*! \brief Random number seed error.
   *
   * Library members throw this exception when a random number
   * generator is provided with an illegitimate starting seed value.
   * For example, the lecuyer class requires seeds within a certain
   * range to operate properly and will throw this exception when
   * seeded with a number outside of that range.
   */
  class scythe_randseed_error:public scythe_exception
  {
  public:
    scythe_randseed_error(const std::string & file,
                      		const std::string & function,
                      		const unsigned int &line,
                      		const std::string & message = "",
                      		const bool & halt = false) throw ()
			: scythe_exception ("SCYTHE RANDOM SEED ERROR", file, function,
					line, message, halt)
    {
    }
  };

  /*! \brief Matrix style error.
   *
   * Library members throw this exception when they are asked to
   * operate on a Matrix of the incorrect style.  Some routines
   * require specifically a concrete Matrix or view to work correctly.
   * For example, only views may reference other matrices; invoking
   * the reference function on a concrete Matrix will generate this
   * exception.
   */
  class scythe_style_error:public scythe_exception
	{
		public:
			scythe_style_error(const std::string& file,
					const std::string& function,
					const unsigned int& line,
					const std::string& message = "",
					const bool& halt = false) throw ()
				:	scythe_exception("SCYTHE STYLE ERROR", file, function,
						line, message, halt)
			{}
	};

  /*! \brief LAPACK Internal Error
   *
   * Library members throw this exception when an underlying LAPACK or
   * BLAS routine indicates that an internal error has occurred.
   * 
   */
  class scythe_lapack_internal_error:public scythe_exception
	{
		public:
			scythe_lapack_internal_error(const std::string& file,
					const std::string& function,
					const unsigned int& line,
					const std::string& message = "",
					const bool& halt = false) throw ()
				:	scythe_exception("SCYTHE LAPACK/BLAS INTERNAL  ERROR", file, 
            function, line, message, halt)
			{}
	};

  /*! \brief Unexpected call to default error.
   *
   * This error should not occur.  If it occurs in your code, please
   * contact the Scythe developers to report the problem.
   * 
   */
  class scythe_unexpected_default_error:public scythe_exception
	{
		public:
			scythe_unexpected_default_error(const std::string& file,
					const std::string& function,
					const unsigned int& line,
					const std::string& message = "",
					const bool& halt = false) throw ()
				:	scythe_exception("SCYTHE UNEXPECTED DEFAULT ERROR", file, 
            function, line, message, halt)
      {}
	};

  // The definition of our terminate handler described above
  inline void scythe_terminate ()
  {
    std::cerr << serr->what() << std::endl;
    std::cerr << std::endl;
    abort ();
  }

}        // end namspace SCYTHE

#endif /* SCYTHE_ERROR_H */
