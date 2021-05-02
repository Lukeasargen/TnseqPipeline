bowtie2: [https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)
---

# Install bowtie2

There are a few ways to get bowtie2. I use the first method.

__Method 1 - apt package manager__

1. Start the Ubuntu app. It will open a terminal.

![open_terminal](/docs/images/open_terminal.png)

2. Run `apt install bowtie`. You will need to enter your password. The terminal spaces will remain empty as you type your password. If you want to also get bowtie2, run `apt install bowtie2`.

![apt_bowtie2](/docs/images/apt_bowtie2.png)

3. Done

__Method 2 - get compiled binaries__

No instructions yet.

## Verify bowtie2 install

The binaries are stored in /bin or /usr/bin and they called:
```
bowtie2
bowtie2-build
bowtie2-inspect
```
Run all three of the commands listed above in any directory of the linux terminal. All commands should display help messages.

Here is the help message from `bowtie2-inspect`:

![message_1](/docs/images/message_1.png)

