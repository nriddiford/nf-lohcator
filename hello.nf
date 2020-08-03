#!/usr/bin/env nextflow

ch = Channel.from(1,2,3)
println(ch)
ch.view()

params.greeting = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

process splitletters {

  input:
  val x from greeting_ch

  output:
  file 'chunk_*' into letters_ch

  """
  printf '$x' | split -b 6 - chunk_
  """
}

process convertToUpper {

  input:
  file y from letters_ch.flatten()

  output:
  stdout into result

  """
  cat $y | tr '[a-z]' '[A-Z]'
  """
}

result.view{ it.trim() }
// result.view()
