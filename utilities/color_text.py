from colorama import init as colorama_init
from colorama import Fore, Style
import sys

colorama_init(strip=False)

def error(text):
  print(Fore.RED + text + Style.RESET_ALL, file=sys.stderr)

def info(text):
  print(Fore.CYAN + text + Style.RESET_ALL, file=sys.stderr)
