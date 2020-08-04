# Groovy exceptions. E.g.:
```{nextflow}
throw new IllegalArgumentException("Unkown arg")
```

# operators on channels
```{nextflow}
.flatten()
.collect() # pulls all channel inputs into list

```
# Docker

`touch Dockerfile`
### add commands in Dockerfile
`docker build -t my-image .`
`docker run my-image`
### -it for interactive mode, where we can 'jump into' the container

### mount filesytem inside container (so you can access files)
`docker run --volume $HOME:$HOME --workdir $PWD -my-image`

### see images

`docker images`

### tag

`docker tag my-image nriviera/my-image:v1`

### push

`docker push nriviera/my-image:v1`

### build
docker build --tag nriviera/bwa-nf-test:test1.0 .

docker push nriviera/bwa-nf-test:test1.0



# resume run from cached results
-resume

# reporting

-with-docker -with-report -with-trace -with-timeline -with-dag dag.png

The -with-report option enables the creation of the workflow execution report. Open the file report.html with a browser to see the report created with the above command.

The -with-trace option enables the create of a tab separated file containing runtime information for each executed task. Check the content of the file trace.txt for an example.

The -with-timeline option enables the creation of the workflow timeline report showing how processes where executed along time. This may be useful to identify most time consuming tasks and bottlenecks. See an example at this link.

Finally the -with-dag option enables to rendering of the workflow execution direct acyclic graph representation. Note: this feature requires the installation of Graphviz in your computer. See here for details.

# run from github
nriddiford/nf-bwa-test
