
# USEFUL FUNCTIONS


# CHECK IF GENELISTS_WHOLEBODY DIRECTORYS CONTAIN THE SAME FILES

genelist_directory1 = "GeneLists-WholeBody"
genelist_directory2 = "GeneLists-WholeBody_new"

compare_GeneLists = function(genelist_directory1, genelist_directory2){
  files1 = tolower(list.files(genelist_directory1))
  files2 = tolower(list.files(genelist_directory2))
  total1 = length(files1)
  total2 = length(files2)
  extra_in_1 = 0
  missing_in_1 = 0
  extrafiles1 = c()
  missingfiles1 = c()
  for (f in files1) {
    if (!f  %in% files2) {
      extra_in_1 = extra_in_1 + 1
      extrafiles1 = append(extrafiles1, f)
    } }

  for (f in files2) {
    if (!f %in% files1) {
      missing_in_1 = missing_in_1 + 1
      missingfiles1 = append(missingfiles1, f)
    }
  }
  names_missing_count = c(paste(genelist_directory1, "total"), paste(genelist_directory2, "total"), paste(genelist_directory1, "extra files"), paste(genelist_directory1, "missing files"))

  missing_counts = c(total1, total2, extra_in_1, missing_in_1)
  names(missing_counts) = names_missing_count
  return (list(missing_counts, extrafiles1, missingfiles1))
}

missingFiles = compare_GeneLists(genelist_directory1, genelist_directory2 )


