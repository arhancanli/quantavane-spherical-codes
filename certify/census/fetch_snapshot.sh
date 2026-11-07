#!/bin/zsh
# Downloads the DSpace@MIT archive of Cohn's table (handle 1721.1/153543, data as of 31 July 2026) politely:
# one request at a time, a pause between files, long back-off on HTTP 429; verifies each file's MD5 against
# the checksum DSpace lists. Usage: fetch_snapshot.sh OUTDIR
out=${1:?outdir}
mkdir -p $out
UA="Mozilla/5.0 (academic data census)"
base="https://dspace.mit.edu/server/api/core/bitstreams"
list=(
"table.txt ce27ecfe-fc3b-4a76-bc83-cf3be69d520e 3adb4968fe9b7b9b21f7be541bcaf270"
"data-1.txt 3ff79008-4daf-4afb-b4ba-83d031a528f9 e12292a6b9343b36cd17f99b6ea6185f"
"data-2.txt 99b97d89-dc8b-4871-9ec0-81a18be48d55 820a60d829ea6c6ee044e9717ad9d191"
"data-3.txt ed04c59e-ae0c-47d8-b938-c4d66b8c6b61 557faa850dfe2d62aba3c0e691caeb1c"
"data-4.txt 56b29399-23a3-42fe-b84c-6cdc305bb2e6 a3f50e4f34d2a2df58031018d97e9874"
"data-5.txt 63f25dfd-e073-417b-94f8-29c5e82cde94 defab0fdbf1c7030f5a427c4da7c250a"
"data-6.txt 82a0e34a-ecf7-424b-a49c-f85408c7972d f0a7938851fbcb9860717adf77067b91"
"data-7.txt e1ebd20e-70a5-4b42-9137-38ae6645a1bb db797dbc943645e036f1582ff85c8b85"
"data-8.txt efdfff23-bca2-4fbf-beba-3d42a241a43e 4e6ad1f2e8fd8fa2273c030ce5930b96"
"data-9.txt 984df0de-59b2-4c17-a8fe-060f30eb3e34 52806c39fe75ab5795185af7b89f91d5"
"data-10.txt cbcf9aa5-d561-476b-aca7-ee166294cbee 81e42121fbb6de24f123d17363370df9"
"data-11.txt a87d4887-f5f5-40b3-949a-72a99c7b9ab6 be355c1f653232a8959c51e3bffdbcc0"
"data-12.txt 9cb5d481-f07d-42d7-91fc-1b91daaf9df5 9cf7c1993a307e4d86fa61d3fc1498fa"
)
sleep 120
for row in $list; do
  parts=(${=row}); name=$parts[1]; id=$parts[2]; md=$parts[3]
  if [[ -f $out/$name ]] && [[ $(md5 -q $out/$name) == $md ]]; then echo "$(date -u +%FT%TZ) have $name"; continue; fi
  for attempt in 1 2 3 4 5 6; do
    code=$(curl -sSL -A "$UA" -o $out/$name.part -w "%{http_code}" "$base/$id/content")
    if [[ $code == 200 ]] && [[ $(md5 -q $out/$name.part) == $md ]]; then
      mv $out/$name.part $out/$name; echo "$(date -u +%FT%TZ) ok $name $(stat -f %z $out/$name) bytes md5 $md"; break
    fi
    echo "$(date -u +%FT%TZ) $name attempt $attempt: http $code, backing off"; rm -f $out/$name.part; sleep $((300*attempt))
  done
  sleep 90
done
echo "$(date -u +%FT%TZ) done"
