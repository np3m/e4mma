#include <yaml-cpp/yaml.h>
#include <yaml-cpp/emitter.h>

#include <sstream>
#include <fstream>

using namespace YAML;

int main(int argc, char *argv[]) {
  YAML::Emitter emitter;
  emitter << "Hello world!";

  std::ofstream fout("file.yaml");
  fout << emitter.c_str();
  return 0;   
}