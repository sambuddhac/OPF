#!/usr/bin/env ruby
require 'benchmark'

# family of different problem sizes
problem_size = [30, 100, 300, 1000, 3000, 10000, 30000,100000,300000]
num_runs = 5

f = File.open("opsp_prof.txt", "w+")
f.puts "# #{(1..num_runs).to_a.join(" ")}"

problem_size.each { |num_nets|  
  f.print "#{num_nets} "
  
  num_runs.times {
    system("./create_network -t 64 -n #{num_nets} -d 0.11")
    result = ""
    measurement = Benchmark.realtime do
      result = `./d_opf -r 1 -t 64 -n 2000`
    end
    
    puts result
    last_line = result.scan(/(\d+) ([\-e0-9\.]+) ([\-e0-9\.]+) ([\-e0-9\.]+) ([\-e0-9\.]+) ([\-e0-9\.]+)/).last
    f.print "#{last_line[0]}(#{measurement}) "
  }
  f.puts
}

f.close()

