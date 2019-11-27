#!/usr/bin/env ruby
require 'rubygems'

require 'net/http'
require 'json'
require 'progressbar'
require 'pp'

Host, Port = 'cvxgen.stanford.edu', 80
SecretKey = 'foobar'

def request command, dict
  begin
    req = Net::HTTP::Post.new("/api/#{command}", {'Content-Type' =>'application/json'})
    dict['key'] = SecretKey
    req.body = dict.to_json

    response = Net::HTTP.start(Host, Port) do |http|
      http.read_timeout = 500 # longer timeouts to wait for codegen
      http.request(req)
    end

    JSON.load response.body
  rescue SocketError
    puts "Error: Failed to connect to #{Host}:#{Port}."
    exit 1
  end
end

def red text
  "\e[31m#{text}\e[0m"
end

$verbose = false

def check prob_desc
  response = request(:check, {"prob_desc" => prob_desc})

  errors = response['errors']
  error_lines = errors.map(&:first)

  prob_desc.lines.each_with_index do |l,i|
    if error_lines.include?(i + 1)
      print red " #{'%3d' % (i+1)}  #{l.rstrip}"
    else
      if $verbose
        puts " #{'%3d' % (i+1)}  #{l}"
      end
    end
  end

  puts
  if errors.empty?
    puts "No errors."
  else
    puts "Errors:"
    errors.each do |i,l|
      puts red " #{'%3d' % i}  #{l}"
    end
    puts
  end

  puts "#{response['kkt_size']}x#{response['kkt_size']} KKT matrix with #{response['nonzeros']} non-zero elements."
end

def generate prob_desc
  puts "starting code generation..."

  response = request(:generate, {"prob_desc" => prob_desc})

  version = response['version']
  puts "Generating problem \##{version}..."
  version
end

def update version
  request(:update, {"version" => version})
end



unless ARGV[0] and File.exists?(ARGV[0])
  puts "Error: Please specify a valid filename."
  exit 1
end

puts "reading #{ARGV[0]}..."
open(ARGV[0]) do |f|
  prob_desc = f.read
  check prob_desc
  version = generate prob_desc

  pending = true
  while pending
    response = update version
    sleep 0.1 # small delay to allow process to start.
    pending = response['pending']
    if (prop = response['permutation_proportion'].to_f) > 0
      begin
        # this should work, but sometimes throws a stack error
        pbar ||= ProgressBar.new 'permutation', 1.0
        pbar.set prop
      rescue Exception => e
        # check if e is systemstackerror
        
        # do nothing if we get a stack error
        # next
      end
    end
    sleep 0.5 if pending
  end
  puts
  case response['state']
  when 'complete'
    puts "Code generation complete."
    puts "Generation time: %.2f s, including permutation time: %.2f s." % [response['codegen_elapsed'], response['perm_elapsed']]
    print "Fill-in: %.2f. " % response['fill_in']
    puts "Code size: %d kB." % response['code_size_in_kb']
    puts "Downloading..."
    Net::HTTP.start(Host, Port) do |http|
      response = http.get "/api/download/#{version}"
      open('/tmp/cvxgen.zip', 'wb') do |f|
        f.write response.body
      end
    end
    puts "File download complete."
    puts "Extracting to cvxgen/"
    `unzip /tmp/cvxgen.zip -d cvxgen/`
    `mv cvxgen/cvxgen/* cvxgen/; rm -r cvxgen/cvxgen`
  else
    puts "Code generation failed."
    pp response
  end
end
