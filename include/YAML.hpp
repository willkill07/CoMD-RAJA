#ifndef YAML_OUTPUT_HPP_
#define YAML_OUTPUT_HPP_

#include <stdio.h>

extern FILE* yamlFile;

void
yamlBegin(void);

void
yamlEnd(void);

void
yamlAppInfo(FILE* file);

void
printSeparator(FILE* file);

#endif
