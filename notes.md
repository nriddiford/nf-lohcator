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

docker build --tag nriviera/bwa-nf-test:test1.0 .

docker push nriviera/bwa-nf-test:test1.0
