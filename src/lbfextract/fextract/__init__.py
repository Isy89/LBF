"""
Fextract module
===============

In the fextract module, the package implements default coverage-based features.

These are the default coverage features implemented:

    - ***fragment coverage***. It takes into account all position of a fragment going from its left most to its right most position.
    - ***midpoint coverage***. It takes into account only the central position of the fragment.
    - ***middle-n points coverage***. It takes into account n positions from the left and right of the midpoint of a fragment.
    - ***coverage around dyads***. It uses only inferred positions of the fragment, at which the dyads may be located
    - ***sliding coverage***. It uses all positions of each fragment from its left most position to its right most position.
    - ***Peter-Ulz (central-60bp) coverage***. It uses only the positions [53, 113) and [-113, -53) positions of reads in a pair 
    - ***windowed protection score (WPS)***. Number of positions derived by fragments completely overlapping a window around a given position. In lbfextract we furhter normalized it based on coverage at each position

More about this topic can be found in the article: `LBFextract: unveiling transcription factor dynamics from liquid biopsy data`,
which can be found `here <link_to_article>`.
"""