#ifndef _FILEIO_H_
#define _FILEIO_H_

#define UNUSED(x) ((void)(x))
#include <cstdio>   //fopen
#include <stdlib.h>  //malloc
#include <string.h>  //strchr
#include <stdarg.h>  // va_list

#include <string>
#include <vector>

// Windows platform does not have good autoconf tools, so skip loading config.h
// _WIN32 is on for both Win32 and Win64
// refer: "On all Windows builds (even 64-bit ones) ‘_WIN32’ will be defined: on 64-bit Windows builds also ‘_WIN64’, and on OS X ‘__APPLE__’ is defined.18"
#ifndef _WIN32
#include "config.h"
#endif

#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

// cannot forward declare an typdef anonymous struct
// http://stackoverflow.com/questions/804894/forward-declaration-of-a-typedef-in-c
// so include the header file
#include "bgzf.h"
#include "R.h"

// #define IO_DEBUG
typedef enum FileType {
  PLAIN = 0,
  GZIP = 1,
#ifdef HAVE_BZIP2
  BZIP2 = 2,
#endif
  BGZIP = 3,
  IO_FILE_TYPE_UNKNOWN = 99
} FileType;

/**
 * Sample usage:
 * std::string line;
 * LineReader* f = AbstractFileReader::open("abc");
 * while (f->readLine(line) > 0) {
 *   ...
 * }
 * delete f;
 */
class AbstractFileReader {
 public:
  typedef enum FileType {
    PLAIN = 0,
    GZIP = 1,
#ifdef HAVE_BZIP2
    BZIP2 = 2,
#endif
#ifdef _USE_KNETFILE
    KNET = 3,
#endif
    IO_FILE_TYPE_UNKNOWN = 99
  } FileType;
  virtual ~AbstractFileReader(){};  // make it virtual so subclass types can
                                    // close file handle
  static AbstractFileReader* open(const char* fileName);
  static void close(AbstractFileReader** f);
  // virtual functions
  // each specific file type will need to implement the following function
  // virtual int readLine(std::string* line) = 0;
  // virtual int readLineBySep(std::vector<std::string>* fields, const char*
  // sep) = 0;
  virtual int getc() = 0;
  virtual bool isEof() = 0;
  virtual void close() = 0;
  virtual int read(void* buf, int len) = 0;
  // common utility function
  static FileType checkFileType(const char* fileName);

 protected:
  AbstractFileReader(){};  // forbid explicit create AbstractFileReader class.
};

class PlainFileReader : public AbstractFileReader {
 public:
  PlainFileReader(const char* fileName) : fp(NULL) {
    this->open(fileName);
#ifdef IO_DEBUG
    REprintf("PlainFileReader() open %s\n", fileName);
#endif
  };
  virtual ~PlainFileReader() {
#ifdef IO_DEBUG
    REprintf("~PlainFileReader() close\n");
#endif
    this->close();
  };

  // get a char, if EOF, return EOF
  int getc() { return ::getc(this->fp); }
  // check eof
  bool isEof() { return (feof(this->fp) != 0); }
  // open
  FILE* open(const char* fileName) {
    this->fp = fopen(fileName, "r");
    if (!this->fp) {
      REprintf("ERROR: Cannot open %s\n", fileName);
    }
    return this->fp;
  }
  // close
  void close() {
    if (this->fp) {
      fclose(fp);
      fp = NULL;
    }
  }
  int read(void* buf, int len) {
    return ::fread(buf, sizeof(char), len, this->fp);
  };

 private:
  FILE* fp;
};

//////////////////////////////////////////////////////////////////////
// Gzip reading class
#include <zlib.h>
class GzipFileReader : public AbstractFileReader {
 public:
  GzipFileReader(const char* fileName) : fp(NULL) {
    this->open(fileName);
#ifdef IO_DEBUG
    REprintf("GzipFileReader() open %s\n", fileName);
#endif
  };
  virtual ~GzipFileReader() {
#ifdef IO_DEBUG
    REprintf("~PlainFileReader() close\n");
#endif
    this->close();
  };

  // get a char, if EOF, return EOF
  int getc() { return gzgetc(this->fp); }
  // check eof
  bool isEof() { return (gzeof(this->fp) != 0); }
  // open
  gzFile open(const char* fileName) {
    this->fp = gzopen(fileName, "r");
    if (!this->fp) {
      REprintf("ERROR: Cannot open %s\n", fileName);
    }
    return this->fp;
  }
  // close
  void close() {
    if (this->fp) {
      gzclose(fp);
      fp = NULL;
    }
  }
  int read(void* buf, int len) { return gzread(this->fp, buf, len); };

 private:
  gzFile fp;
};

//////////////////////////////////////////////////////////////////////
// Bzip2 reading class
#ifdef HAVE_BZIP2

#include <bzlib.h>
class Bzip2FileReader : public AbstractFileReader {
 public:
  Bzip2FileReader(const char* fileName) : fp(NULL) {
    this->open(fileName);
#ifdef IO_DEBUG
    REprintf("Bzip2FileReader() open %s\n", fileName);
#endif
  };
  virtual ~Bzip2FileReader() {
#ifdef IO_DEBUG
    REprintf("~Bzip2FileReader() close\n");
#endif
    if (this->fp) {
      BZ2_bzclose(fp);
    }
  };

  // get a char, if EOF, return EOF
  int getc() {
    char c;
    this->bzerror = BZ_OK;
    int nBuf = BZ2_bzRead(&this->bzerror, this->bzp, &c, sizeof(char));
    UNUSED(nBuf);
    if (this->bzerror == BZ_OK) {
      return c;
    } else {
      return EOF;
    }
  }
  // check eof
  bool isEof() { return (this->bzerror == BZ_STREAM_END); }
  // open
  BZFILE* open(const char* fileName) {
    this->fp = fopen(fileName, "rb");
    if (!this->fp) {
      REprintf("ERROR: Cannot open %s\n", fileName);
      return NULL;
    }
    this->bzp = BZ2_bzReadOpen(&this->bzerror, this->fp, 0, 0, NULL, 0);

    if (this->bzerror != BZ_OK) {
      BZ2_bzReadClose(&bzerror, this->bzp);
      REprintf("ERROR: Cannot open %s\n", fileName);
      return NULL;
    }
    return this->bzp;
  }
  // close
  void close() {
    if (this->bzerror != BZ_STREAM_END) {
      BZ2_bzReadClose(&this->bzerror, this->bzp);
      /* omit code to handle error */
    } else {
      BZ2_bzReadClose(&this->bzerror, this->bzp);
    }
    if (this->fp) fclose(this->fp);
    this->fp = NULL;
    this->bzp = NULL;
    this->bzerror = 0;
  };
  int read(void* buf, int len) {
    return BZ2_bzRead(&this->bzerror, this->bzp, buf, len);
  };

 private:
  FILE* fp;
  BZFILE* bzp;
  int bzerror;
};
#endif

#ifdef _USE_KNETFILE
// need to include this before R.h to avoid redefinition of Realloc
// #include "knetfile.h"
class KnetFileReader : public AbstractFileReader {
 public:
  KnetFileReader(const char* fileName) : fp(NULL) {
#ifdef _WIN32
    knet_win32_init();
#endif
    this->open(fileName);
#ifdef IO_DEBUG
    REprintf("KnetFileReader() open %s\n", fileName);
#endif
  };
  virtual ~KnetFileReader() {
#ifdef IO_DEBUG
    REprintf("~KnetFileReader() close\n");
#endif
    this->close();
#ifdef _WIN32
    knet_win32_destroy();
#endif
  };

  // get a char, if EOF, return EOF
  int getc() {
    // return bgzf_getc(this->fp);
    int c;
    if (knet_read(this->fp, &c, sizeof(char)) == 0) {
      return EOF;
    } else {
      return c;
    }
  }
  // check eof
  bool isEof() {
    // we don't know the exact file size over internet;
    // and we only know eof from read() return value
    if (!this->fp) { // probably files are not openned
      return true;
    }
    return (this->eof);
    // return bgzf_check_EOF(this->fp);
  }
  // open
  knetFile* open(const char* fileName) {
    this->eof = false;
    this->fp  = knet_open(fileName, "r");
    //     fpr = knet_open(fileName, "rb");
// #ifdef _WIN32    
//     fpr = knet_open(fileName, "rb");
// #else
//     fpr = knet_open(fileName, "rb");
// #endif
//     if (!fpr) {
//       REprintf("ERROR: Cannot open %s\n", fileName);
//     }

    // this->fp = bgzf_read_init();
    // this->fp->fp = fpr;

    // this->fp = bgzf_open(fileName, "r");
    if (!this->fp) {
      REprintf("ERROR: Cannot open %s\n", fileName);
    }
    return this->fp;
  }
  // close
  void close() {
    if (this->fp) {
      // bgzf_close(this->fp);
      knet_close(this->fp);
      fp = NULL;
    }
    this->eof = true;
  }
  int read(void* buf, int len) {
    // return bgzf_read(this->fp, buf, len);
    off_t l = knet_read(this->fp, buf, len);
    if (l < len || l == 0) {
      this->eof =  true;
    }
    return l;
  }

 private:
  // BGZF* fp;
  knetFile* fp;
  bool eof;  // when read() fewer bytes or zero
};
#endif
//////////////////////////////////////////////////////////////////////

/**
 * Example code:

 BufferedReader br("Makefile", 200);
 char buf[500] = {};
 int nRead = 0;
 while ( (nRead = br.read(buf, 500)) > 0) {
 for (int i = 0; i < nRead; i++) {
 printf("%c", buf[i]);
 }
 }
*/
class BufferedReader : public AbstractFileReader {
 public:
  BufferedReader(const char* fileName, int bufferCapacity)
      : bufCap(0), bufEnd(0), bufPtr(0), buf(NULL), fp(NULL) {
#ifdef IO_DEBUG
    REprintf("BufferedReader open %s\n", fileName);
#endif
    // initialize buf
    if (bufferCapacity == 0) {
      REprintf(
          "Buffer size should be greater than 0, now use default buffer size "
          "8096 instead of %d.\n",
          bufferCapacity);
      this->bufCap = 8096;
    } else {
      this->bufCap = (int)(bufferCapacity);
    }
    this->buf = new char[this->bufCap];
    if (!this->buf) {
      REprintf("Cannot allocate buffer for BufferedReader.\n");
      return;
    }
    this->bufPtr = 0;
    this->bufEnd = 0;
    // initialize fp
    this->fp = AbstractFileReader::open(fileName);
    if (!this->fp) {
      REprintf("Canont open file %s\n", fileName);
      this->fp = NULL;
    }
  }
  virtual ~BufferedReader() { this->close(); }
  int getc() {
    if (this->bufPtr == this->bufEnd) {  // buffer all used, need to refresh
      this->bufEnd = this->fp->read(this->buf, this->bufCap);
      this->bufPtr = 0;
    }

    if (this->bufPtr < this->bufEnd)
      return (this->buf[this->bufPtr++]);
    else
      return EOF;
  }
  bool isEof() {
    if (this->fp && this->fp->isEof() && this->bufPtr == this->bufEnd) {
      return true;
    }
    return false;
  }
  void close() {
#ifdef IO_DEBUG
    REprintf("BufferedReader close\n");
#endif
    // delete fp
    if (this->fp) {
      AbstractFileReader::close(&fp);
    }
    this->fp = NULL;
    // delete buf
    if (this->buf) {
      delete[] this->buf;
      this->buf = NULL;
      this->bufCap = 0;
      this->bufPtr = 0;
      this->bufEnd = 0;
    }
    this->buf = NULL;
  }
  int read(void* buf_, int len) {
    // use current buffer to fill in buf
    int idx = 0;
    while (this->bufPtr < this->bufEnd && len > 0) {
      ((char*)buf_)[idx++] = this->buf[this->bufPtr++];
      len--;
    }
    if (len == 0) {
      return idx;
    }
    // fill rest of buf
    int nRead = this->fp->read(((char*)buf_) + idx, len);
    idx += nRead;
    // refill buffer
    this->bufEnd = this->fp->read(this->buf, this->bufCap);
    this->bufPtr = 0;
    return idx;
  }

 private:
  int bufCap;  // capacity of the buffer
  int bufEnd;  // bufPtr should not read beyond bufEnd(incluive)
  int bufPtr;  // from which buffer begins to read
  char* buf;
  AbstractFileReader* fp;
};

/** Example code:
 // LineReader lr(fn);
 // while(lr.readLine(&line)>0){
 //         Rprintf( "%s\n", line.c_str());
 // }
 */
class LineReader {
 public:
  LineReader(const std::string& fileName) { init(fileName.c_str()); }
  LineReader(const char* fileName) { init(fileName); }
  void init(const char* fileName) {
#ifdef IO_DEBUG
    REprintf("LineReader open %s\n", fileName);
#endif
    this->fp = new BufferedReader(fileName, 1024);
    if (!this->fp) {
      REprintf("Canont open file %s\n", fileName);
      this->fp = NULL;
    }
  }
  LineReader(AbstractFileReader* fp_) { this->fp = fp_; }
  virtual ~LineReader() {
    if (this->fp) {
      fp->close();
      delete this->fp;
      this->fp = NULL;
    }
#ifdef IO_DEBUG
    REprintf("LineReader close\n");
#endif
  }
  // return number of characters read.
  // when reading an empty line, will return 1, as we read '\n', however, line
  // will be empty
  // when reading the end, we will return 0
  int readLine(std::string* line) {
    // //assert(this->fp && line);
    if (this->fp->isEof()) return 0;
    line->clear();
    char c;
    unsigned nRead = 0;
    while (true) {
      c = this->fp->getc();
      if (c == EOF) {
        return nRead;
      } else if (c == '\r') {
        // skip this
        continue;
      } else if (c == '\n') {
        ++nRead;
        return nRead;
      } else {  // normal characters
        ++nRead;
        line->push_back(c);
      }
    }
    // //assert(false); // should not reach here
    return 0;
  };
  // return number of fields read.
  // when reading an empty line, will return 1, meaning 1 field are read,
  // although its content is empty
  // consecutive separators, e.g. \t\t, will yield empty field
  // when reading to the EOF, will return 0.
  int readLineBySep(std::vector<std::string>* fields, const char* sep) {
    // //assert(this->fp && fields && sep);
    if (this->fp->isEof()) return 0;
    fields->clear();
    char c;
    std::string s;
    while (true) {
      c = this->fp->getc();
      if (c == EOF) {
        fields->push_back(s);
        return fields->size();
      } else if (c == '\r') {
        // skip this
        continue;
      } else if (c == '\n') {
        fields->push_back(s);
        return fields->size();
      } else if (strchr(sep, c) != NULL) {  // separator
        fields->push_back(s);
        s.clear();
      } else {  // normal characters
        s.push_back(c);
      }
    }
    // //assert(false); // should not reach here
    return 0;
  };

 private:
  AbstractFileReader* fp;
};

//////////////////////////////////////////////////////////////////////
// FileWriter related classes
class AbstractFileWriter {
 public:
  /// when open is successful, return 0; else: return non-zero
  virtual int open(const char* fn, bool append = false) = 0;
  virtual void close() = 0;
  virtual int write(const char* s) = 0;
  virtual int writeLine(const char* s) = 0;
  virtual ~AbstractFileWriter() = 0;
};

class TextFileWriter : public AbstractFileWriter {
 public:
  TextFileWriter(const char* fn, bool append = false) {
    if (this->open(fn, append)) {
      REprintf("Cannot create text file %s\n", fn);
    }
  }
  virtual ~TextFileWriter() {
#ifdef IO_DEBUG
    REprintf("TextFileWriter desc()\n");
#endif
    this->close();
  }
  int open(const char* fn, bool append = false) {
    if (append)
      this->fp = fopen(fn, "a");
    else
      this->fp = fopen(fn, "w");
    if (!this->fp) {
      REprintf("ERROR: Cannot open %s for write\n", fn);
      return -1;
    }
    return 0;
  }
  void close() {
    if (this->fp) {
      fclose(this->fp);
      this->fp = NULL;
    }
  };
  int write(const char* s) { return fputs(s, this->fp); };
  int writeLine(const char* s) {
    int ret = fputs(s, this->fp);
    fputc('\n', this->fp);
    return (ret + 1);
  };
  int printf(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    int ret = vfprintf(fp, fmt, args);
    va_end(args);
    return ret;
  };

 private:
  FILE* fp;
};  // end TextFileWriter

class GzipFileWriter : public AbstractFileWriter {
 public:
  GzipFileWriter(const char* fn, bool append = false) {
    if (this->open(fn, append)) {
      REprintf("Cannot create gzip file %s\n", fn);
    }
  }
  virtual ~GzipFileWriter() {
    this->close();
#ifdef IO_DEBUG
    REprintf("GzipFileWriter desc()\n");
#endif
  };
  int open(const char* fn, bool append = false) {
    if (append) REprintf("Gzip does not support appending.\n");
    this->fp = gzopen(fn, "wb");
    if (!this->fp) {
      REprintf("ERROR: Cannot open %s for write\n", fn);
      return -1;
    }
    return 0;
  }
  void close() {
    if (this->fp) {
      gzclose(this->fp);
      this->fp = NULL;
    }
  };
  int write(const char* s) { return gzputs(this->fp, s); };
  int writeLine(const char* s) {
    int ret = gzputs(this->fp, s);
    gzputc(this->fp, '\n');
    return (ret + 1);
  };

 private:
  gzFile fp;
};  // end GzipFileWriter

#ifdef HAVE_BZIP2
class Bzip2FileWriter : public AbstractFileWriter {
 public:
  Bzip2FileWriter(const char* fn, bool append = false) {
    if (this->open(fn, append)) {
      REprintf("Cannot create bzip2 file %s\n", fn);
    }
  }
  virtual ~Bzip2FileWriter() {
    this->close();
#ifdef IO_DEBUG
    REprintf("Bzip2FileWriter desc()\n");
#endif
  };
  int open(const char* fn, bool append = false) {
    if (append) REprintf("bzip2 does not support appending.\n");
    this->fp = fopen(fn, "wb");
    if (fp == NULL) return -1;

    this->bzp = BZ2_bzWriteOpen(
        &this->bzerror, this->fp, 9, 0,
        30);  // block size is 9, 0 means silent, 30 means working factor
    if (this->bzerror != BZ_OK) {
      BZ2_bzWriteClose(&bzerror, this->bzp, 0, 0, 0);  // 0: abandon, 0: results
                                                       // of # of bytes for
                                                       // input, 0: results of #
                                                       // of bytes outputted.
      REprintf("ERROR: Cannot open %s for write\n", fn);
      return -1;
    }
    return 0;
  }
  void close() {
    BZ2_bzWriteClose(&bzerror, this->bzp, 0, 0, 0);
    if (bzerror != BZ_OK) {
    };
    if (this->fp) fclose(this->fp);

    this->bzp = NULL;
    this->fp = NULL;
  };
  int write(const char* s) {
    int ret = strlen(s);
    BZ2_bzWrite(&this->bzerror, this->bzp, (void*)s, ret);
    if (this->bzerror != BZ_OK) {
      this->close();
      return -1;
    }
    return ret;
  };
  int writeLine(const char* s) {
    int ret = strlen(s);
    BZ2_bzWrite(&this->bzerror, this->bzp, (void*)s, ret);
    if (this->bzerror != BZ_OK) {
      this->close();
      return -1;
    }
    char buf[] = "\n";
    BZ2_bzWrite(&this->bzerror, this->bzp, buf, 1);
    if (this->bzerror != BZ_OK) {
      this->close();
      return -1;
    }
    return (ret + 1);
  };

 private:
  FILE* fp;
  BZFILE* bzp;
  int bzerror;
};  // end Bzip2FileWriter
#endif

class BGZipFileWriter : public AbstractFileWriter {
 public:
  BGZipFileWriter(const char* fn, bool append = false) {
    if (this->open(fn)) {
      REprintf("Cannot create BGzip file %s\n", fn);
    }
  }
  virtual ~BGZipFileWriter() {
    this->close();
#ifdef IO_DEBUG
    REprintf("BGZipFileWriter desc()\n");
#endif
  };
  /**
   * @param append: ignored
   */
  int open(const char* fn, bool append = false);
  void close();
  int write(const char* s);
  int writeLine(const char* s);

 private:
  BGZF* fp;
};  // end BGZipFileWriter

#define DEFAULT_WRITER_BUFFER 4096
class BufferedFileWriter : public AbstractFileWriter {
 public:
  BufferedFileWriter(AbstractFileWriter* f_,
                     int bufLen_ = DEFAULT_WRITER_BUFFER) {
    this->bufLen = DEFAULT_WRITER_BUFFER;
    this->buf = new char[bufLen_ + 1];  // last char in the buffer is always '\0'
    // that help to use fputs()
    if (!this->buf) {
      REprintf("Cannot create BufferedFileWriter\n");
      REprintf("Critical error happening!\n");  // abort();
    }
    this->buf[bufLen_] = '\0';
    this->bufPtr = 0;

    if (!this->buf) {
      REprintf("Buffer allocation failed!\n");
    }
    this->f = f_;
  }
  ~BufferedFileWriter() {
    if (this->buf) {
      delete[] this->buf;
      this->buf = NULL;
    }
#ifdef IO_DEBUG
    REprintf("BufferedFileWriter desc()\n");
#endif
  };
  int open(const char* fn, bool append = false) {
    return this->f->open(fn, append);
  };
  void close() {
    this->flush();
    // this->f->close();
  };
  int write(const char* s) {
    int nbyte = 0;
    int i = 0;
    while (s[i] != '\0') {
      this->buf[this->bufPtr++] = s[i++];
      nbyte++;
      if (this->bufPtr == this->bufLen) {
        this->f->write(this->buf);
        this->bufPtr = 0;
      }
    }
    return nbyte;
  };
  int writeLine(const char* s) {
    int ret = this->write(s);
    this->write("\n");
    return (ret + 1);
  };
  int flush() {
    this->buf[this->bufPtr] = '\0';
    this->f->write(this->buf);
    this->bufPtr = 0;
    return 0;
  };

 private:
  char* buf;
  int bufLen;
  int bufPtr;
  AbstractFileWriter* f;
};  // end BufferedFileWriter

/**
 * design:
 *  a high level file class, underlying using BufferedFileWriter
 * usage:
 * FileWriter* fout = new FileWriter("a.txt", "w");
 * fout->write("abc");
 * fout->writeLn("abc");
 * fout->close();
 * delete fout->write;
 */
class FileWriter {
 public:
  FileWriter(const char* fileName, bool append = false) {
    // int l = strlen(fileName);
    if (this->checkSuffix(fileName, ".gz")) {
      this->fpRaw = new GzipFileWriter(fileName, append);
#ifdef HAVE_BZIP2
    } else if (this->checkSuffix(fileName, ".bz2")) {
      this->fpRaw = new Bzip2FileWriter(fileName, append);
#endif
    } else {
      this->fpRaw = new TextFileWriter(fileName, append);
    }
    this->fp = new BufferedFileWriter(this->fpRaw);
    if (!this->fpRaw || !this->fp) {
      REprintf("Cannot create file\n");
      REprintf("Critical error happening!\n");  // abort();
    }

    this->createBuffer();
  }
  FileWriter(const char* fileName, FileType t) {
    bool append = false;
    if (PLAIN == t) {
      this->fpRaw = new TextFileWriter(fileName, append);
    } else if (GZIP == t) {
      this->fpRaw = new GzipFileWriter(fileName, append);
#ifdef HAVE_BZIP2
    } else if (BZIP2 == t) {
      this->fpRaw = new Bzip2FileWriter(fileName, append);
#endif
    } else if (BGZIP == t) {
      this->fpRaw = new BGZipFileWriter(fileName, append);
    } else {
      REprintf("Unrecognized file type, use plain text format instead!\n");
      this->fpRaw = new TextFileWriter(fileName, append);
    }

    this->fp = new BufferedFileWriter(this->fpRaw);
    if (!this->fpRaw || !this->fp) {
      REprintf("Cannot create file\n");
      REprintf("Critical error happening!\n");  // abort();
    }

    this->createBuffer();
  };
  void createBuffer() {
    // create buffer for formatted string
    this->bufLen = 1024;
    this->buf = new char[this->bufLen];
    if (!this->buf) {
      REprintf("Cannot allocate printf buffer for FileWriter.\n");
    };
  };
  void close() {
    if (this->fp) {
      this->fp->close();
      delete this->fp;
      this->fp = NULL;
    }
    if (this->fpRaw) {
      delete this->fpRaw;
      this->fpRaw = NULL;
    }
    if (this->buf) {
      delete[] this->buf;
      this->buf = NULL;
    }
#ifdef IO_DEBUG
    REprintf("FileWriter desc()\n");
#endif
  };
  ~FileWriter() { this->close(); };
  int write(const char* s) { return this->fp->write(s); };
  int write(const std::string& s) { return this->fp->write(s.c_str()); }
  int writeLine(const char* s) {
    int ret = this->fp->write(s);
    this->fp->write("\n");
    return (ret + 1);
  };
  // if @param fileName ends with @param suffix, then return true;
  static bool checkSuffix(const char* fileName, const char* suffix) {
    int lf = strlen(fileName);
    int ls = strlen(suffix);
    if (lf < ls) return false;
    for (int i = lf - ls, j = 0; j < ls;) {
      if (fileName[i++] != suffix[j++]) return false;
    }
    return true;
  };

  /**
   * format string to this->buf, then write it out
   */
  int printf(const char* fmt, ...) {
    // we'll put the formatted string  to internal buffer
    va_list args;
    int ret;
    int newBufLen;

    while (1) {
      /* Try to print in the allocated space. */
      va_start(args, fmt);
      ret = vsnprintf(this->buf, this->bufLen, fmt, args);
      va_end(args);
      /* If that worked, return the string. */
      if (ret > -1 && ret < this->bufLen) {
        return this->write(this->buf);
      }
      /* Else try again with more space. */
      if (ret > -1)             /* glibc 2.1 */
        newBufLen = ret + 1;    /* precisely what is needed */
      else                      /* glibc 2.0 */
        newBufLen = bufLen * 2; /* twice the old size */
      increaseBufferTo(newBufLen);
    }
  };

  /*   while (( ret = vsnprintf(this->buf, this->bufLen, fmt, args)) >=
   * this->bufLen){ */
  /*     this->increaseBufferTo( ret + this->bufLen + 1); */
  /*   }; */
  /*   va_end(args); */
  /*   return this->write(this->buf); */
  /* }; */

 private:
  void increaseBufferTo(int newBufLen) {
    delete[] this->buf;
    this->buf = new char[newBufLen];
    if (!this->buf) {
      REprintf("Cannot increase printf buffer for FileWriter.\n");
      REprintf("Critical error happening!\n");  // abort();
    };
    this->bufLen = newBufLen;
  };
  AbstractFileWriter* fp;
  AbstractFileWriter* fpRaw;
  char* buf;
  int bufLen;
};  // end class FileWriter

#endif /* _FILEIO_H_ */
