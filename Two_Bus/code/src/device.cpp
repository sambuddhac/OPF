#include "device.hpp"

Device::Device(int N, int len) : len(len), N(N), terminals(new Terminal*[N]())
{
  // DC terminals first
  for(int i = 0; i < N; i++)
    terminals[i] = new Terminal(len);
}

Device::~Device()
{
  if(terminals) { 
    for(int i = 0; i < N; i++)
      if(terminals[i]) delete terminals[i];
    delete[] terminals;
  }
}

void Device::print()
{
  for(int i = 0; i < N; i++) {
    std::cout << terminals[i] << " ";
    std::cout << *terminals[i] << std::endl;
  }
}