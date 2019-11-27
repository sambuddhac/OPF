#!/usr/bin/env ruby

def insert_namespace filename, namespace
  tmp = File.open("tmp", 'w')
  File.open(filename, 'r') do |f|
    while(line = f.gets)
      if line.match(/^\/\*|^\n|^#(.*)+/)
        tmp.puts line
      else 
        break
      end
    end
    tmp.puts "namespace #{namespace} {"
    tmp.puts line
    while(line = f.gets)

      tmp.puts line
    end
    tmp.puts "}"
  end
  tmp.close
  # replace anything matching "#endif\n}" with "}\n#endif"
  text = File.read("tmp")
  File.open("tmp", "w") { |f|
    f.puts text.gsub(/#endif\n\}/, "}\n#endif")
  }
  
  `cp tmp #{filename}`
  `rm tmp`
  
  
end

def insert_using filename, namespace
  tmp = File.open("tmp", 'w')
  File.open(filename, 'r') do |f|
    while(line = f.gets)
      if line.match(/^\/\*|^\n|^#(.*)+/)
        tmp.puts line
      else 
        break
      end
    end
    tmp.puts "using namespace #{namespace};"
    tmp.puts line
    while(line = f.gets)
      tmp.puts line
    end
  end
  tmp.close
  `cp tmp #{filename}`
  `rm tmp`
end

has_namespace = true
unless ARGV[0]
  puts "No namespace encapsulation will be used."
  has_namespace = false
end

if has_namespace
  directory = ARGV[1]
  unless ARGV[1]
    puts "Looking under default cvxgen directory."
    directory = "cvxgen"
  end

  unless File.exists?(directory)
    puts "Error: please specify a valid directory."
    exit 1
  end

  puts "Inserting namespace '#{ARGV[0]}'"
  Dir[directory + "/*.c"].each do |f|
    if f.include? "testsolver"
      insert_using f, ARGV[0]

      # scope the function load_default_data in the namespace
      # also scope the variables
      text = File.read(f)
      File.open(f, "w") { |solver|
        text = text.gsub(/^Vars vars;/, "Vars #{ARGV[0]}::vars;");
        text = text.gsub(/^Params params;/, "Params #{ARGV[0]}::params;");
        text = text.gsub(/^Workspace work;/, "Workspace #{ARGV[0]}::work;");
        text = text.gsub(/^Settings settings;/, "Settings #{ARGV[0]}::settings;");
        
        solver.puts text.gsub(/void load_default_data/, "namespace #{ARGV[0]} {\nvoid load_default_data")
        solver.puts "}"
      }
    else
      insert_namespace f, ARGV[0]
    end
    # everybody else inserts "namespace ARGV[1] { ... }"
    # testsolver inserts "using namespace ARGV[1]"
  end
  
  # modify makefile to use g++
  text = File.read("cvxgen/Makefile")
  File.open("cvxgen/Makefile", "w") { |f|
    f.puts text.gsub(/gcc/, "g++")
  }
  
  # scope objects in to proper namespace
  insert_namespace "cvxgen/solver.h", ARGV[0]
  
  # rename the #defines for the header file to allow multiple inclusions 
  # of solver.h
  text = File.read("cvxgen/solver.h")
  File.open("cvxgen/solver.h", "w") { |f|
    text = text.gsub(/^#ifndef SOLVER_H/, "#ifndef " + 
      "SOLVER_#{ARGV[0]}_H".upcase);
    f.puts text.gsub(/^#define SOLVER_H/, "#define " + 
      "SOLVER_#{ARGV[0]}_H".upcase);
  }
  
end