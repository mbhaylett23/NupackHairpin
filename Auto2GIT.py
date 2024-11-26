#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:03:31 2024

@author: michaelbh23
"""

import subprocess
import os

# Correct WSL path
path = "/mnt/c/Users/michael.bruyns/NupackGIT/Nupack4GIT/MainCode"

# Change the directory
try:
    os.chdir(path)
    print(f"Current working directory: {os.getcwd()}")
except FileNotFoundError:
    print(f"Error: The directory '{path}' does not exist.")
except PermissionError:
    print(f"Error: Permission denied to access '{path}'.")

# Function to run a shell command
def run_command(command):
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")

# Configure Git user.name and user.email
user_name = "Michael Bruyns-Haylett"  # Replace with your name
user_email = "mbruynshaylett@.gmail.com"  # Replace with your email

run_command(["git", "config", "--global", "user.name", user_name])
run_command(["git", "config", "--global", "user.email", user_email])



# Personal Access Token (PAT)
pat = "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"  # Replace with your actual PAT

# Update the remote URL with the PAT
remote_url = f"https://{pat}@github.com/mbhaylett23/NupackHairpin.git"
run_command(["git", "remote", "set-url", "origin", remote_url])

# Push changes
run_command(["git", "push", "origin", "main"])
# Add all changes
run_command(["git", "add", "."])

# Commit changes with a message
commit_message = input("Enter commit message: ")
# commit_message = "Auto-update from Spyder"
run_command(["git", "commit", "-m", commit_message])

# Push changes to the remote repository
run_command(["git", "push", "origin", "main"])