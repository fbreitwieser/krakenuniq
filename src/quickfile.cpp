/*
 * Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "quickfile.hpp"

using std::string;

namespace kraken {

QuickFile::QuickFile() {
  valid = false;
  fptr = NULL;
  filesize = 0;
  fd = -1;
}

QuickFile::QuickFile(string filename_str, string mode, size_t size) {
  open_file(filename_str, mode, size);
}

void QuickFile::open_file(string filename_str, string mode, size_t size, bool lock) {
  const char *filename = filename_str.c_str();
  int o_flags = mode == "w"
                  ? O_RDWR | O_CREAT | O_TRUNC
                  : mode == "r" ? O_RDONLY : O_RDWR;
  int m_flags = mode == "r" ? MAP_PRIVATE : MAP_SHARED;

  fd = open(filename, o_flags, 0666);
  // Second try for R/W if failure was due to non-existence
  if (fd < 0 && mode == "rw" && errno == ENOENT) {
    o_flags |= O_CREAT;
    fd = open(filename, o_flags, 0666);
  }
  if (fd < 0)
    err(EX_OSERR, "unable to open %s", filename);

  if (o_flags & O_CREAT) {
    if (lseek(fd, size - 1, SEEK_SET) < 0)
      err(EX_OSERR, "unable to lseek (%s)", filename);
    if (write(fd, "", 1) < 0)
      err(EX_OSERR, "write error (%s)", filename);
    filesize = size;
  }
  else {
    struct stat sb;
    if (fstat(fd, &sb) < 0)
      err(EX_OSERR, "unable to fstat %s", filename);
    filesize = sb.st_size;
  }

  fptr = (char *)mmap(0, filesize, PROT_READ | PROT_WRITE, m_flags, fd, 0);
  if (fptr == MAP_FAILED)
    err(EX_OSERR, "unable to mmap %s", filename);

  if (lock) {
    int ret = mlock(fptr, filesize);
    if (ret) {
      if (errno == ENOMEM) {
        err(EX_OSERR, "unable to memory lock %s due to soft limits. Please set ulimit -l higher", filename);
      }
      else {
        err(EX_OSERR, "unable to memory lock %s. Is there enough available memory?", filename);
      }
    }
  }
  valid = true;

  if (mode == "r") {
    size_t page_size = sysconf(_SC_PAGESIZE);
    size_t n_pages = (filesize + page_size - 1) / page_size;

    unsigned char* buffer = new unsigned char[n_pages];
    mincore(fptr, filesize, buffer);
    size_t cached_pages = 0;
    for (size_t i = 0; i < n_pages; ++i) {
      if (buffer[i] & 1) {
        cached_pages++;
      }
    }
    double cached_gb = (double) cached_pages * page_size / pow(1024, 3);
    double filesize_gb = (double) n_pages * page_size / pow(1024, 3);
    std::cerr << "Database file " << filename_str << " - " << cached_pages << "/" << n_pages << " of pages cached, "
              << std::fixed << std::setprecision(2) << cached_gb << "/" << filesize_gb << " GiB, "
              << (float) cached_pages / n_pages * 100 << "%" << std::endl;
  }
}

void QuickFile::load_file() {
  int thread_ct = 1;
  int thread = 0;
  #ifdef _OPENMP
  int old_thread_ct = omp_get_max_threads();
  if (old_thread_ct > 4)
    omp_set_num_threads(4);
  thread_ct = omp_get_max_threads();
  #endif

  size_t page_size = getpagesize();
  char buf[thread_ct][page_size];

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

#ifdef _OPENMP
    #pragma omp for schedule(dynamic)
#endif
    for (size_t pos = 0; pos < filesize; pos += page_size) {
      size_t this_page_size = filesize - pos;
      if (this_page_size > page_size)
        this_page_size = page_size;
      memcpy(buf[thread], fptr + pos, this_page_size);
    }
  }

  #ifdef _OPENMP
  omp_set_num_threads(old_thread_ct);
  #endif
}

char * QuickFile::ptr() {
  return valid ? fptr : NULL;
}

size_t QuickFile::size() {
  return valid ? filesize : 0;
}

QuickFile::~QuickFile() {
  close_file();
}

void QuickFile::sync_file() {
  msync(fptr, filesize, MS_SYNC);
}

void QuickFile::close_file() {
  if (! valid)
    return;
  sync_file();
  munmap(fptr, filesize);
  close(fd);
  valid = false;
}

// from http://programanddesign.com/cpp/human-readable-file-size-in-c/
char* readable_fs(double size/*in bytes*/, char *buf) {
    int i = 0;
    const char* units[] = {"B", "kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
    while (size > 1024) {
        size /= 1024;
        i++;
    }
    sprintf(buf, "%.*f %s", i, size, units[i]);
    return buf;
}



std::vector<char> slurp_file(std::string filename, size_t lSize) {
  FILE * pFile;
  size_t result;

  pFile = fopen ( filename.c_str() , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  if (lSize > 0) {
    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);
  }
  
  char buf[50];
  readable_fs(lSize, buf);
  std::cerr << "Getting " << filename << " into memory (" << buf << ") ..."; 

  // copy the file into the vector:
  std::vector<char> buffer(lSize);
  result = fread (buffer.data(),1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  fclose (pFile);

  std::cerr << " Done" << std::endl;
  return(buffer);
}



} // namespace
