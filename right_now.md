
# To Do

## Work through Pemba lineages example

- [ ] Test each function

- [ ] Profile and see if there are any speed ups/simplifications

## Check parallelization

- [ ] Improvements locally

- [ ] Improvements cluster (i.e. up the number of cores?)

- [ ] Still reproducible?

- [x] explicitly managing DT threads [actually okay because upon fork, i.e. foreach
      becomes single threaded and swithces back after fork is closed!]

## Other bits
- [x] Add extraction of date & time difference between cases to consensus links
- [ ] Make creating graphs more flexible/easier
- [ ] Interactive + labeled graphs with tooltips (vizgraph, linked to map + ts if possible!) 
- [ ] Automake gifs (map + timeseries)
