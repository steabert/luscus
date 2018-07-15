/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<unistd.h>
#include<sys/types.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

MOL *mol;
FILE *backup = NULL;
int backup_count = 0;

void close_backup(void)
{
  if (backup) fclose(backup);
  backup = NULL;
}

void remove_backup(void)
{
  char *backup_filename = get_backup_filename();
  backup_count = 0;
  unlink(backup_filename);
  free(backup_filename);
}

void open_backup(int app_backup)
{
  char *backup_filename;
  backup_filename = get_backup_filename();

  if (app_backup) backup = fopen(backup_filename, "a");
  else backup = fopen(backup_filename, "r+");

  free(backup_filename);
}

char *get_backup_filename(void)
{
  char *input_filename;
  char *extension;
  char *backup_filename;
  char *path, *filename, *filename_without_extension;
 
  /*deconstruct path*/
  input_filename = strdup(get_input_filename());
  extension = get_extension_pointer(input_filename);
  path = get_path_without_file(input_filename);
  filename = get_filename_pointer_from_path(input_filename);
  filename_without_extension = get_filename_without_extension(filename);

  /*allocate memory*/
  backup_filename = (char*) malloc(sizeof(char) * (strlen(path) + strlen(filename_without_extension) + strlen(input_filetypes[0].extension) + 3));

  /*construct new path*/
  backup_filename[0] = 0;
  strcat(backup_filename, path);
  strcat(backup_filename, ".");
  strcat(backup_filename, filename_without_extension);
  strcat(backup_filename, ".");
  strcat(backup_filename, input_filetypes[0].extension);

  /*free memory*/
  free(input_filename);
  free(path);
  free(filename_without_extension);

  return backup_filename;
}

void append_backup(void)
{
  if (backup == NULL) open_backup(1);
  write_backup();
  close_backup();
  backup_count++;
}

void write_backup(void)
{
  write_all_data_to_output(backup, 0);
  fprintf(backup, " <BREAK>\n");
}

char *get_extension_pointer(char *filename)
{
  char *tmp;
  tmp = strrchr(filename, 46);
  if (tmp) return tmp;
  else return filename + strlen(filename); /*no extension, return pointer to the end*/
}

char *get_filename_pointer_from_path(char *path)
{
  char *tmp;
#ifdef WINDOWS
  tmp = strrchr(path, 92);
#else
  tmp = strrchr(path, 47);
#endif
  if (tmp) return tmp+1;
  else return path;
}

char *get_path_without_file(char *path0)
{
  int n, i;
  char *path1 = NULL;
  char *file_pointer;
  char null = 0;
#ifdef WINDOWS
  file_pointer = strrchr(path0, 92);
#else
  file_pointer = strrchr(path0, 47);
#endif
  if (!file_pointer) return strdup(&null);

  path1 = (char*) malloc(sizeof(char) * (file_pointer - path0 + 2));
  strncpy(path1, path0, file_pointer - path0 + 1);
  path1[file_pointer - path0 + 1] = 0;

  return path1;
}

char *get_filename_without_extension(char *filename)
{
  char *filename_without_extension;
  char *ext = get_extension_pointer(filename);

  filename_without_extension = malloc(sizeof(char) * (ext - filename + 1));
  strncpy(filename_without_extension, filename, ext - filename);
  filename_without_extension[ext - filename] = 0;
  return filename_without_extension;
}

void get_backup(void)
{
  int i;
  off_t len;
  char *line;
  int next = 1;
  int *m_ngrids = NULL;
  if (backup_count <= 1)
  {
    luscus_gtk_push_message_to_statusbar2("Allready the oldest change!");
    return;
  }

  m_ngrids = (int*) malloc(n_geometries * sizeof(int));
  for(i = 0; i < n_geometries; i++)
    m_ngrids[i] = mol[i].ngrids;

  if (backup == NULL) open_backup(0);
  backup_search_last_backup_point();
/*  delete_all_orbitals_from_list();*/

  parse_gv_file(backup);

  read_all_grids_from_all_sections(m_ngrids);
  /*--------------*/
  luscus_gtk_show_or_hide_widgets();
  /*search for next "break"*/
  while(next && !feof(backup))
  {
    line = read_line(backup);
    if (strcasestr(line, "<break>") != 0) next = 0;
    free(line);
  }

  len = ftell(backup);
  if(ftruncate(fileno(backup), len)) fprintf(stderr, "ERROR: can't truncate backup file\n");

  close_backup();
  backup_count--;
  free(m_ngrids);
}

void backup_search_last_backup_point(void)
{
  int ibackup = 0;
  char *line;

  rewind(backup);
  while(ibackup < backup_count - 2)
  {
    line = read_line(backup);
    if (strcasestr(line, "<BREAK>")) ibackup++;
    free(line);
    if (feof(backup))
    {
      luscus_gtk_push_message_to_statusbar2("Error in reading backup data!");
      break;
    }
  }
}

