#!/bin/bash

echo "[-] Removing all tc rules from lo"

tc qdisc del dev lo root 2>/dev/null

tc qdisc show dev lo

