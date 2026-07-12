#!/usr/bin/env python3
from __future__ import annotations

import argparse

from common import deep_update, read_yaml, write_yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--defaults", required=True)
    parser.add_argument("--dataset-config", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    defaults = read_yaml(args.defaults)
    dataset = read_yaml(args.dataset_config) if args.dataset_config else {}
    write_yaml(args.output, deep_update(defaults, dataset))


if __name__ == "__main__":
    main()
