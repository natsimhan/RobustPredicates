# RobustPredicates

PHP port of [RobustPredicates][robustpredicates], a fast robust predicates for
computational geometry.

Provides reliable 2D ~~and 3D~~ point orientation tests (orient2d, ~~orient3d,
incircle, insphere~~) that are not susceptible to floating point errors (without
sacrificing performance).

:warning: Only orient2d is ported to PHP for the moment.


## Installation

The preferred method of installation is via [Composer][]. Run the following
command to install the package and add it as a requirement to your project's
`composer.json`:

```bash
composer require natsimhan/robust-predicates
```


## Copyright and License

This code was placed in the public domain by its original author,
[Jonathan Richard Shewchuk][jonathanshewchuk]. You may use it as you see fit,
but attribution is appreciated. Please see [LICENSE][] for more information.


[robustpredicates]: https://github.com/mourner/robust-predicates
[composer]: http://getcomposer.org/
[jonathanshewchuk]: https://people.eecs.berkeley.edu/~jrs/
[license]: ./LICENSE
