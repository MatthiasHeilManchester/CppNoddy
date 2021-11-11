// hacky include; stolen from oomph-lib; no include guards and stuff. 

#include <map>

  //====================================================================
  // Namespace for command line arguments
  //====================================================================
  namespace CommandLineArgs
  {
    /// Structure to store information on a command line argument
    template<class T>
    struct ArgInfo
    {
      /// Proper constructor
      ArgInfo(const bool& is_set, T* arg_pt, const std::string& doc)
      {
        this->is_set = is_set;
        this->arg_pt = arg_pt;
        this->doc = doc;
      }

      /// Default constructor. We need this to be able to store these things in
      /// maps.
      ArgInfo()
      {
        this->is_set = false;
        this->arg_pt = 0;
        this->doc = "";
      }

      // /// Assignment operator. We need this to be able to store these things
      // in
      // /// maps.
      // void operator=(const ArgInfo& that)
      //  {
      //   this->is_set = that.is_set;
      //   this->arg_pt = that.arg_pt;
      //   this->doc = that.doc;
      //  }

      /// Has this argument been set?
      bool is_set;

      /// The place to put the argument's value when it is set
      T* arg_pt;

      /// Information about what the argument does
      std::string doc;
    };

    // Number of arguments + 1
    extern int Argc;

    // Arguments themselves
    extern char** Argv;

    /// Map to indicate an input flag as having been set
    extern std::map<std::string, ArgInfo<bool>> Specified_command_line_flag;

    /// Map to associate an input flag with a double -- specified via pointer
    extern std::map<std::string, ArgInfo<double>>
      Specified_command_line_double_pt;

    /// Map to associate an input flag with an int -- specified via pointer
    extern std::map<std::string, ArgInfo<int>> Specified_command_line_int_pt;

    /// Map to associate an input flag with an unsigned -- specified via pointer
    extern std::map<std::string, ArgInfo<unsigned>>
      Specified_command_line_unsigned_pt;

    /// Map to associate an input flag with a string -- specified via pointer
    extern std::map<std::string, ArgInfo<std::string>>
      Specified_command_line_string_pt;

    // Set values
    extern void setup(int argc, char** argv);

    // Doc the command line arguments
    extern void output();

    /// Specify possible argument-free command line flag
    extern void specify_command_line_flag(
      const std::string& command_line_flag,
      const std::string& documentation = "Undocumented");

    /// Specify possible command line flag that specifies a double,
    /// accessed via pointer
    extern void specify_command_line_flag(
      const std::string& command_line_flag,
      double* arg_pt,
      const std::string& documentation = "Undocumented");

    /// Specify possible command line flag that specifies an int,
    /// accessed via pointer
    extern void specify_command_line_flag(
      const std::string& command_line_flag,
      int* arg_pt,
      const std::string& documentation = "Undocumented");

    /// Specify possible command line flag that specifies an unsigned,
    /// accessed via pointer
    extern void specify_command_line_flag(
      const std::string& command_line_flag,
      unsigned* arg_pt,
      const std::string& documentation = "Undocumented");

    /// Specify possible command line flag that specifies a string,
    /// accessed via pointer
    extern void specify_command_line_flag(
      const std::string& command_line_flag,
      std::string* arg_pt,
      const std::string& documentation = "Undocumented");

    /// Check if specified flag has been set (the associated
    /// value will have been assigned directly)
    extern bool command_line_flag_has_been_set(const std::string& flag);

    /// Document the values of all flags (specified or not)
    extern void doc_all_flags(
     std::ostream& outstream); // hierher = *std::cout.stream_pt());

    /// Document specified command line flags
    extern void doc_specified_flags();

    /// Document available command line flags
    extern void doc_available_flags();

    /// Helper function to check if command line index is legal
    extern void check_arg_index(const int& argc, const int& arg_index);

    /// Parse command line, check for recognised flags and assign
    /// associated values
    extern void parse_and_assign(
      int argc, char* argv[], const bool& throw_on_unrecognised_args = false);

    /// Parse previously specified command line, check for
    /// recognised flags and assign associated values
    extern void parse_and_assign(
      const bool& throw_on_unrecognised_args = false);

  } // namespace CommandLineArgs


  //====================================================================
  /// Namespace for command line arguments
  //====================================================================
  namespace CommandLineArgs
  {
    /// Number of arguments + 1
    int Argc;

    /// Arguments themselves
    char** Argv;

    /// Map to indicate an input flag as having been set
    std::map<std::string, ArgInfo<bool>> Specified_command_line_flag;

    /// Map to associate an input flag with a double -- specified via pointer
    std::map<std::string, ArgInfo<double>> Specified_command_line_double_pt;

    /// Map to associate an input flag with an int -- specified via pointer
    std::map<std::string, ArgInfo<int>> Specified_command_line_int_pt;

    /// Map to associate an input flag with an unsigned -- specified via pointer
    std::map<std::string, ArgInfo<unsigned>> Specified_command_line_unsigned_pt;

    /// Map to associate an input flag with a string -- specified via pointer
    std::map<std::string, ArgInfo<std::string>>
      Specified_command_line_string_pt;

    /// Set values
    void setup(int argc, char** argv)
    {
      Argc = argc;
      Argv = argv;
    }

    /// Doc the command line arguments
    void output()
    {
      std::cout << "You are running the program: " << CommandLineArgs::Argv[0]
                 << std::endl;
      std::cout << "with the following command line args: " << std::endl;
      std::stringstream str;
      for (int i = 1; i < CommandLineArgs::Argc; i++)
      {
        str << CommandLineArgs::Argv[i] << " ";
      }
      std::cout << str.str() << std::endl;
    }


    /// Specify possible argument-free command line flag
    void specify_command_line_flag(const std::string& command_line_flag,
                                   const std::string& doc)
    {
      Specified_command_line_flag[command_line_flag] =
        ArgInfo<bool>(false, 0, doc);
    }

    /// Specify possible command line flag that specifies a double,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   double* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_double_pt[command_line_flag] =
        ArgInfo<double>(false, arg_pt, doc);
    }

    /// Specify possible command line flag that specifies an int,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   int* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_int_pt[command_line_flag] =
        ArgInfo<int>(false, arg_pt, doc);
    }

    /// Specify possible command line flag that specifies an unsigned,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   unsigned* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_unsigned_pt[command_line_flag] =
        ArgInfo<unsigned>(false, arg_pt, doc);
    }

    /// Specify possible command line flag that specifies a string,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   std::string* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_string_pt[command_line_flag] =
        ArgInfo<std::string>(false, arg_pt, doc);
    }


    /// Check if command line flag has been set (value will have been
    /// assigned directly).
    bool command_line_flag_has_been_set(const std::string& flag)
    {
      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        if (it->first == flag)
        {
          return it->second.is_set;
        }
      }

      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      return false;
    }

    /// Document the values of all flags (specified or not).
    void doc_all_flags(std::ostream& outstream)
    {
      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        outstream << it->first << " " << it->second.is_set << std::endl;
      }
      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        outstream << it->first << " " << *(it->second.arg_pt) << std::endl;
      }
      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        outstream << it->first << " " << *(it->second.arg_pt) << std::endl;
      }
      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        outstream << it->first << " " << *(it->second.arg_pt) << std::endl;
      }
      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        // Quote blank strings, otherwise trying to parse the output in any way
        // will go wrong.
        std::string arg_string = *(it->second.arg_pt);
        if (arg_string == "")
        {
          arg_string = "\"\"";
        }

        outstream << it->first << " " << arg_string << std::endl;
      }
    }

    /// Document specified command line flags
    void doc_specified_flags()
    {
      std::cout << std::endl;
      std::cout << "Specified (and recognised) command line flags:\n";
      std::cout << "----------------------------------------------\n";

      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        if (it->second.is_set)
        {
          std::cout << it->first << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          std::cout << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          std::cout << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          std::cout << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          std::cout << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      std::cout << std::endl;
    }


    /// Document available command line flags
    void doc_available_flags()
    {
      std::cout << std::endl;
      std::cout << "Available command line flags:\n";
      std::cout << "-----------------------------\n";

      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        std::cout << it->first << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        std::cout << it->first << " <double> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        std::cout << it->first << " <int> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        std::cout << it->first << " <unsigned> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        std::cout << it->first << " <string> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      std::cout << std::endl;
    }


    /// Helper function to check if command line index is legal
    void check_arg_index(const int& argc, const int& arg_index)
    {
      if (arg_index >= argc)
      {
        output();
        doc_available_flags();
        std::stringstream error_stream;
        error_stream
          << "Tried to read more command line arguments than\n"
          << "specified. This tends to happen if a required argument\n"
          << "to a command line flag was omitted, e.g. by running \n\n"
          << "     ./a.out -some_double \n\n rather than\n\n"
          << "     ./a.out -some_double 1.23 \n\n"
          << "To aid the debugging I've output the available\n"
          << "command line arguments above.\n";
        std::cout << error_stream.str() << std::endl;
        abort();
      }
    }


    /// Parse command line, check for recognised flags and assign
    /// associated values
    void parse_and_assign(int argc,
                          char* argv[],
                          const bool& throw_on_unrecognised_args)
    {
      // Keep looping over command line arguments
      int arg_index = 1;
      while (arg_index < argc)
      {
        bool found_match = false;

        if ((strcmp(argv[arg_index], "-help") == 0) ||
            (strcmp(argv[arg_index], "--help") == 0))
        {
          std::cout
            << "NOTE: You entered --help or -help on the command line\n";
          std::cout << "so I'm going to tell you about this code's\n";
          std::cout << "available command line flags and then return.\n";
          doc_available_flags();

#ifdef OOMPH_HAS_MPI
          int flag;
          MPI_Initialized(&flag);
          if (flag != 0) MPI_Helpers::finalize();
#endif
          std::cout<< "Shutting down...\n";
          exit(0);
        }


        // Check if the flag has been previously specified as a simple argument
        // free command line argument
        for (std::map<std::string, ArgInfo<bool>>::iterator it =
               Specified_command_line_flag.begin();
             it != Specified_command_line_flag.end();
             it++)
        {
          if (it->first == argv[arg_index])
          {
            Specified_command_line_flag[argv[arg_index]].is_set = true;
            found_match = true;
            break;
          }
        }

        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) a double
          for (std::map<std::string, ArgInfo<double>>::iterator it =
                 Specified_command_line_double_pt.begin();
               it != Specified_command_line_double_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the double. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = atof(argv[arg_index]);
              found_match = true;
              break;
            }
          }
        }


        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) an int
          for (std::map<std::string, ArgInfo<int>>::iterator it =
                 Specified_command_line_int_pt.begin();
               it != Specified_command_line_int_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the int. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = atoi(argv[arg_index]);
              found_match = true;
              break;
            }
          }
        }


        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) an
          // unsigned
          for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
                 Specified_command_line_unsigned_pt.begin();
               it != Specified_command_line_unsigned_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the unsigned. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = unsigned(atoi(argv[arg_index]));
              found_match = true;
              break;
            }
          }
        }


        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) a string
          for (std::map<std::string, ArgInfo<std::string>>::iterator it =
                 Specified_command_line_string_pt.begin();
               it != Specified_command_line_string_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the string. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = argv[arg_index];
              found_match = true;
              break;
            }
          }
        }


        // Oh dear, we still haven't found the argument in the list.
        // Maybe it was specified wrongly -- issue warning.
        if (!found_match)
        {
          // Construct the error message
          std::string error_message = "Command line argument\n\n";
          error_message += argv[arg_index];
          error_message += "\n\nwas not recognised. This may be legal\n";
          error_message += "but seems sufficiently suspicious to flag up.\n";
          error_message += "If it should have been recognised, make sure you\n";
          error_message += "used the right number of \"-\" signs...\n";

          if (throw_on_unrecognised_args)
          {
            error_message += "Throwing an error because you requested it with";
            error_message += " throw_on_unrecognised_args option.";
            cout << error_message << std::endl;
            // throw OomphLibError(
            //   error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
            abort();
          }
          else
          {
           // Create an Oomph Lib warning
           cout << error_message << std::endl;
           // OomphLibWarning(
           //   error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
           abort();
          }
        }

        arg_index++;
      }
    }


    /// Parse previously specified command line, check for
    /// recognised flags and assign associated values
    void parse_and_assign(const bool& throw_on_unrecognised_args)
    {
      parse_and_assign(Argc, Argv, throw_on_unrecognised_args);
    }

  } // namespace CommandLineArgs

