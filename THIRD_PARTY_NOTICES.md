# Third-party attribution

The three-contact reduction uses the facet-selection lemma of Andrii Arman, Andriy Bondarenko, Andriy Prymak and Danylo Radchenko, *Borsuk number in four dimensions is at most seven*, arXiv:2610.39078 (2026).

The interval arithmetic in `src/contact/six_balls.py` is derived from their `seven_balls_single.py`, available at <https://github.com/andriyprm/4d7balls>. The contact parameterization and interval enclosure in `src/contact/hd_enclose.py` use this arithmetic. Attribution is also retained in the source headers and the article.

NumPy, SciPy (including Qhull and HiGHS), and pycddlib/cddlib are external dependencies and are not bundled in this repository. They retain their respective authorship and licenses. No additional rights to third-party materials are granted by this notice.
