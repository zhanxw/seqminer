#include "MmapFile.h"

#include <cstdint>

size_t getFileSize(const char* fileName) {
  struct stat fileStat;
  int fd = open(fileName, O_RDONLY);
  if (fd == -1) {
    REprintf("Cannot open file");
    // exit(1);
    return SIZE_MAX;
  }

  if (fstat(fd, &fileStat)) {
    REprintf("Cannot fstat() file");
    // exit(1);
    return SIZE_MAX;
  }
  close(fd);
  return fileStat.st_size;
}

int MmapFile::open(const char* fileName) {
  this->filedes = ::open(fileName, O_RDONLY);
  if (filedes < 0) {
    REprintf("Cannot open file");
    // exit(1);
    return -1;
  }
  this->fileSize = ::getFileSize(fileName);
  if (data) {
    this->close();
  }
  this->data = mmap(0, this->fileSize, PROT_READ, MAP_SHARED, this->filedes, 0);
  if (this->data == MAP_FAILED) {
    REprintf("mmap() failed!");
    // exit(1);
    return -1;
  }
  return 0;
}

void MmapFile::close() {
  munmap(this->data, this->fileSize);
  this->data = 0;
}
