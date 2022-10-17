# Overview

> Performed analysis for the Group Project concerning the first year of PhD in Data Science funded by 
[SFI Centre for Research Training in Foundations of Data Science](https://www.data-science.ie/).

# Authors

- [Andrê F. B. Menezes](https://andrmenezes.github.io/about/)
- [Beatrice Tropea](https://www.data-science.ie/user/beatrice+tropea/)
- [Cormac Monaghan](https://www.data-science.ie/user/cormac+monaghan/)


# **Goals for the Project**
(Group Meeting - 23/09/2022)


#### **Goal 1**
Conducting multiple comparisons between groups using a robust or more complex statistical model\
Using an empirical Bayes to estimate the variance in the denominator of the T-Test\
Potential ways of doing this include:

* Utilizing the Limma package in bioconductor
* Utilizing the package [msqrob2](http://bioconductor.org/packages/release/bioc/html/msqrob2.html) 
  
  
#### **Goal 2**
Identify the highly variable proteins (main proteins) based off their expression\
Potential for networking
  
#### **Goal 3**
Estimating the protein network based off their correlations/interactions

* We think the simplest way to do this would be to use the **Pearson Linear Correlation Coefficient**

Using this we could perform a **cluster analysis**\
Proteins of interest include:

* CspA
* CspE
* DeaD
* LsrF
* YjiM

We could also investigate FliC, along with **persister cell formation pattern**

#### **Other potential goals**
Differences in resistance between ampicillin and ciprofloxacin\
Differences between β-lactamases and ciprofloxacin
  

<div align="center">**Visualization through heat maps**</div>


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
  
  



