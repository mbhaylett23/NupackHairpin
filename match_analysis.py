#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:03:38 2024

@author: michaelbh23
"""

# match_analysis.py

def score_match(subject, query, subject_start, query_start, length):
    '''Compute similarity score between two sequences'''
    score = sum(1 if subject[i + subject_start] == query[i + query_start] else -1 for i in range(length))
    return score

def try_all_matches(subject, query, score_limit):
    '''Find all matches of `subject` in `query` with score >= score_limit'''
    matches = []
    subject_len = len(subject)
    max_start = len(query) - subject_len + 1
    for query_start in range(max_start):
        score = score_match(subject, query, 0, query_start, subject_len)
        if score >= score_limit:
            query_end = query_start + subject_len
            matches.append((True, query_start, query_end, subject_len, score))
    return matches if matches else [(False, None, None, None, None)]