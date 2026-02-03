#input: list of GSM numbers 2 download
#searches entrez database for corresponding SRR numbers and downloads
#concats em if necessary (hopefully)
#  NCBI recommends that users post no more than three URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure to comply with this policy may result in an IP address being blocked from accessing NCBI

base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

for gsm in $*
do
  [[ "$gsm" =~ GSM[0-9]+ ]] || continue
  testy=$(wget -nv -O - "$base""esearch.fcgi?db=sra&RetMax=999&term=""$gsm""&usehistory=yes" 2>/dev/null)

  #have to extract webenv and query_key for future queries
  guid=$(echo $testy | grep -oP "(?<=<Id>)\d+(?=</Id>)" | tac)
  webenv=$(echo $testy | grep -oP "(?<=<WebEnv>).+(?=</WebEnv)")
  querykey=$(echo $testy | grep -oP "(?<=<QueryKey>)\d+(?=</QueryKey>)")

  #j is the GUID we can use to search the SRA database and the SRR numberzzzzz (～﹃～)
  #search for filename="SRR\d+"
  sleep 1
  req=$(wget -nv -O - "$base""efetch.fcgi?db=sra&WebEnv=""$webenv""&Query_key=""$querykey" 2>/dev/null)
  # do IDs one at a time to keep track of which thing to concat
  # Remove duplicates with sort -u
  srrs=$(echo $req | grep -oP '(?<=filename=")SRR\d+' | sort -u)
  echo "For $guid/$gsm:"
  # Download SRRs sequentially (no &) to avoid too many files at once
  for n in $srrs
  do
    echo "  Downloading: $n"
    fasterq-dump -O $gsm $n
    sleep 1
  done

  # Concatenate FASTQ files
  echo "  Concatenating FASTQ files for $gsm..."
  if compgen -G "$gsm/SRR*_1.fastq" > /dev/null; then
    # Paired-end: concatenate _1 and _2 separately
    cat "$gsm"/SRR*_1.fastq > "$gsm"_1.fastq
    cat "$gsm"/SRR*_2.fastq > "$gsm"_2.fastq

    # Compress immediately after concatenation
    echo "  Compressing $gsm paired-end files..."
    gzip "$gsm"_1.fastq &
    gzip "$gsm"_2.fastq &
    wait
  else
    # Single-end: concatenate all
    cat "$gsm"/SRR*.fastq > $gsm.fastq

    # Compress immediately after concatenation
    echo "  Compressing $gsm single-end file..."
    gzip $gsm.fastq
  fi

  # Remove intermediate SRR files to save space
  echo "  Cleaning up intermediate files for $gsm..."
  rm -r $gsm

  # IMPORTANT: ncbi gets grumpy and bans your ip if you send > 3 requests per second
  # have to apply for an api key to fix it
  sleep 1
done

echo ""
echo "============================================"
echo "All downloads and compression complete! ヾ(•ω•\`)o bye"
echo "============================================"
