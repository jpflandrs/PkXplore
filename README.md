# PkXplore


julia --project
using GenieFramework; Genie.loadapp(); up()

export GENIE_ENV=prod
docker build . -t pkxplore 
docker run  -p 8000:8000 pkxplore

