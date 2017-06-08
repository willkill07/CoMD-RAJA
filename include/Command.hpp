#ifndef MYCOMMAND_H
#define MYCOMMAND_H

#include <stdio.h>

struct Command {
  char potDir[1024]{"pots"};
  char potName[1024]{""};
  char potType[1024]{"funcfl"};
  int doeam{0};
  int nx{20};
  int ny{20};
  int nz{20};
  int xproc{1};
  int yproc{1};
  int zproc{1};
  int nSteps{100};
  int printRate{10};
  double dt{1.0};
  double lat{-1.0};
  double temperature{600.0};
  double initialDelta{0.0};

  Command(int argc, char *argv[]);

  Command() = delete;

  Command(Command const &) = default;

  Command(Command &&) = default;

  Command &
  operator=(Command const &) = default;

  void
  printYAML(FILE *);
};

#endif
