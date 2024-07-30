# v0.0 - initial release

# v0.1 - Now with more work!

- doAllQC.py has been reworked to rely on the Entrez API for taxonomy IDs of downloaded sequences, not on parsing the species name from the sequence label.
- We keep sp., cf., etc. species only if we don't have any other species in that genus - previously we removed them completely.
- Added a new COI download script, currently running.
