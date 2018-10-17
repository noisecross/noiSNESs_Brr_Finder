# noiSNESs_Brr_Finder

<p align="center">
<img width="105" src="https://avatars2.githubusercontent.com/u/22390927"></img>
</p>

This tool is strongly based on the good old **SNESSOR** from [Infe](https://www.romhacking.net/community/430) and includes also elements of the amazing **BRRTools** from [Bregalad](https://www.romhacking.net/community/1067).

The **noiSNESs_Brr_Finder** takes a file which contains BRR data, like an SPC or a SNES rom file, and try to rip the data out. The differences with the SNESSOR are:

 - Every possible BRR file is statistically analyzed to discriminate and not to export false positives.
 - The BRR samples can be exported as WAV but also as raw BRR. This makes easy for a romhacker to find where are the sound info chunks in every game and also to swap BRR between games with no losses of quality.
 - The conversion from BRR to WAV is performed using techniques shared by Bregalad like the Gaussian and the treble-boost filters.
