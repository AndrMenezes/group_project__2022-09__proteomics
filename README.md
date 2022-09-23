## Data set

```
Rows: 14,175
Columns: 5
$ protein__id <chr> "ACQ41973.1;P62552", "ACQ41973.1;P62552", "ACQ41973.1;P62552", "ACQ41973.1;P62…
$ group       <chr> "Ampicillin", "Ampicillin", "Ampicillin", "Cefotaxime", "Cefotaxime", "Cefotax…
$ replicate   <int> 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, …
$ variable    <chr> "lqf", "lqf", "lqf", "lqf", "lqf", "lqf", "lqf", "lqf", "lqf", "lqf", "lqf", "…
$ value       <dbl> 24.7033, 23.5168, 26.3269, 25.5626, 26.1038, 26.5440, 25.2792, 26.2555, 25.369…
```


## Ideas to visualize

- Is it interesting to identify highly variable proteins based on their expression?
  - A common approach to selection is based on the protein expression across the population (cultures of bacteria)
  - An usual practice is decompose the protein variability into technical and
  biological by means
  - A plot of mean versus variance could help to visualize the proteins expressions across the samples.
  


## Ideas to analyze

- Estimate and identify the most correlated proteins by group. 

- Use a robust or more complex statistical model to test the most expressed
proteins in each group compared to the control.
  - Package [msqrob2](http://bioconductor.org/packages/release/bioc/html/msqrob2.html) can deal with this.



