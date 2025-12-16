#pragma once
#include "TSystem.h"
#include <iostream>
#include <sstream>

namespace rad{
  namespace files{
    
    // Run shell command and capture output
    std::vector<std::string> runCommand(const std::string& cmd) {
      std::vector<std::string> lines;
      FILE* pipe = gSystem->OpenPipe(cmd.c_str(), "r");
      if (!pipe) {
	std::cerr << "Error: cannot run command " << cmd << std::endl;
        return lines;
      }
      char buffer[1024];
      while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::string line(buffer);
        line.erase(line.find_last_not_of(" \n\r\t") + 1); // trim whitespace
        if (!line.empty()) lines.push_back(line);
      }
      gSystem->ClosePipe(pipe);
      return lines;
    }
    
    // Return list of ROOT files from XRootD directory
    std::vector<std::string> GetXRootDFiles(const std::string redirector,
					    const std::string xrdfsPath,
					    const std::string extension = "edm4eic.root",
					    int maxFiles = -1) {
      std::ostringstream cmd;
      cmd << "xrdfs " << redirector << " ls " << xrdfsPath;
      auto files = runCommand(cmd.str());
    
      // Sort to make sure 0001, 0002, ... order
      std::sort(files.begin(), files.end());
    
      std::vector<std::string> filelist;
      for (auto& f : files) {
        if (TString(f).EndsWith(extension.c_str())) {
	  filelist.push_back("root://" + std::string(redirector) + f);
	  if (maxFiles > 0 && (int)filelist.size() >= maxFiles) break;
	}
      }
    
      if (filelist.empty()) {
        std::cerr << "Warning: No " << extension
                  << " files found under " << xrdfsPath << std::endl;
      }

      return filelist;
    }
    
  }
}
